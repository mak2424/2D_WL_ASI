#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <map>
#include <random>
#include <gsl/gsl_histogram2d.h>

using namespace std;

#define L 5 // linear size of system
#define N L*L // number of spins in the lattice
#define NEIGHCOUNT 4 // number of neighbours of single spin
#define POINTSMAX 100 // maximal (approx) number of possible directions of magnetisation vector
#define EMIN (-N*3-1)
#define EMAX (N*3+1)
#define DOSSIZE1 ulong(EMAX-EMIN) // number of bins in dos histogramm
#define MMAX (N+1)
#define MMIN (-N-1)
#define DOSSIZE2 ulong(MMAX-MMIN) // number of bins in dos histogramm


static double pts[POINTSMAX*3+6]; // array of values of magnetisation vectors
static unsigned points=0; // actual count of possible stated, usually different with NEIGHCOUNT
static unsigned neights[N*NEIGHCOUNT]; // array of neighboring spins
static unsigned parts[N]; // array of particles

static double accuracy=0.8;
static double f=1;
static unsigned walks=22;
static double j=1., d=0.6;
static int rseed;
static unsigned k;

static gsl_histogram2d *g;
static int h[DOSSIZE1*DOSSIZE2];
static char zeros[DOSSIZE1*DOSSIZE2];

static double eNew, eOld, de,
    mNew, mOld, dm,
    mxOld,myOld,mzOld;

inline uint idx(uint a, uint b){ return a*L+b; }

double hamiltonianHeisenbergDMI(
        double ax, double ay, double az,
        double bx, double by, double bz,
        char dirX = 0, char dirY = 0 ) //direction for calculating DMI: basis directions
{
    return -j * (ax*bx+ay*by+az*bz)
            - d * dirX *
            (ay*bz-az*by)
            - d * dirY *
            (az*bx-ax*bz);
}

double hamiltonianHeisenberg(double ax, double ay, double az,
                             double bx, double by, double bz,
                             char dirX = 0, char dirY = 0)
{
    (void)dirX; (void)dirY;
    return -j * (ax*bx+ay*by+az*bz);
}

bool isFlat(){
    unsigned num=0, count=0; double threshold=0, mean=0;
    for (num=0; num<DOSSIZE1*DOSSIZE2; ++num){
        if (zeros[num]) {
            ++count;
            mean += (h[num]-mean)/count;
        }
    }

    threshold = mean * accuracy;

    for (num=0; num<DOSSIZE1*DOSSIZE2; ++num){
        if (zeros[num] && h[num] < threshold){ //критерий плоскости
            //cout<<"# isflat failed at E="<<(num-DOSSHIFT)*PRECISIONFACTOR_E<<"("<<num<<"): "<<h[num]<<" | "<<threshold<<endl;
            return false;
        }
    }
    return true;
}

void norm(){
    unsigned num=0;
    double gmax=gsl_histogram2d_max_val(g);
    double sum=0, a;

    for (num=0; num < DOSSIZE1*DOSSIZE2; ++num){
        if (zeros[num])
            sum+=exp(g->bin[num]-gmax);
    }

    a = gmax + log(sum);

    for (num=0; num < DOSSIZE1*DOSSIZE2; ++num){
        if (zeros[num]){
            g->bin[num] -= a;
        }
    }
}

void dump(){
    ostringstream fname;
    fname<<"dump_"<<N<<"_"<<rseed<<"_"<<k<<".dat";
    ofstream f(fname.str());
    f<<"# WL calculations for 3D spins on flat square lattice"<<endl;
    f<<"# with exchange and DMI interaction"<<endl;
    f<<"# J="<<j<<"\t and D="<<d<<endl;
    f<<"# random seed="<<rseed<<endl;
    f<<"# bins number="<<DOSSIZE1<<" , "<<DOSSIZE2<<endl;
    f<<"# points count="<<POINTSMAX<<"("<<points<<")"<<endl;
    f<<"# emin , emax = "<<EMIN<<" , "<<EMAX<<endl;
    f<<"# mmin , mmax = "<<MMIN<<" , "<<MMAX<<endl;

    for (ulong num1=0; num1<DOSSIZE1; ++num1){
        for (ulong num2=0; num2<DOSSIZE2; ++num2){
            if (zeros[num1*DOSSIZE2+num2]){
                f<<
                    num1<<"\t"<<
                    num2<<"\t"<<
                    g->xrange[num1]<<"\t"<<
                    g->yrange[num2]<<"\t"<<
                    g->bin[num1*DOSSIZE2+num2]<<"\t"<<
                    h[num1*DOSSIZE2+num2]<<endl;
            }
        }
    }
    f.close();
}

void dumpCfg(double e){
    ostringstream fname;
    fname<<"parts_"<<N<<"_"<<rseed<<"_"<<k<<".dat";
    ofstream ff(fname.str());
    ff<<"# config for E="<<e<<endl;
    for (unsigned i=0; i<N; ++i){
        ff<<
             pts[parts[i]*3+0]<<"\t"<<
             pts[parts[i]*3+1]<<"\t"<<
             pts[parts[i]*3+2]<<"\t"<<
             parts[i]<<"\t"<<
             i<<endl;
    }
    ff.close();
    cout<<"# conf dumped"<<endl;
}

void dumpPts(){
    ostringstream fname;
    fname<<"points_"<<N<<"_"<<rseed<<"_"<<k<<".dat";
    ofstream f(fname.str());
    for (unsigned num=0; num<points; ++num){
        f<<pts[num*3+0]<<"\t"<<pts[num*3+1]<<"\t"<<pts[num*3+2]<<endl;
    }
    f.close();
}

void genPoints(){
    //add first point up
    //pts[points*3 + 0] = 0.; //x
    //pts[points*3 + 1] = 0.; //y
    //pts[points*3 + 2] = 1.; //z
    //points+=1;

    double a,d,dTheta,dPhi,theta,phi;
    int mTheta, mPhi, m, n;

    a = 4. * M_PI / POINTSMAX;
    d = sqrt(a);
    mTheta = int(round(M_PI/d));
    dTheta = M_PI/mTheta;
    dPhi = a/dTheta;
    for (m=0; m<mTheta; ++m){
        theta = M_PI * (m+0.5) / mTheta;
        mPhi = int(round(2. * M_PI * sin(theta) / dPhi));

        for (n=0; n<mPhi; ++n){
            phi = 2. * M_PI * n / mPhi;

            pts[points*3 + 0] = sin(theta) * cos(phi); //x
            pts[points*3 + 1] = sin(theta) * sin(phi); //y
            pts[points*3 + 2] = cos(theta); //z

            points+=1;
        }
    }

    // add last point down
    //pts[points*3 + 0] = 0.; //x
    //pts[points*3 + 1] = 0.; //y
    //pts[points*3 + 2] = -1.; //z
    //points+=1;
}

void drop(){
    uint i,j;
    for (i=0;i<L;++i){
        for (j=0;j<L;++j){
            parts[idx(i,j)]=(j%2==0)?0:int(POINTSMAX/2.);

            neights[idx(i,j)*NEIGHCOUNT + 0] = (j==L-1) ? idx(i,0) : idx(i,j+1); //right
            neights[idx(i,j)*NEIGHCOUNT + 1] = (i==0) ? idx(L-1,j) : idx(i-1,j); //up
            neights[idx(i,j)*NEIGHCOUNT + 2] = (j==0) ? idx(i,L-1) : idx(i,j-1); //left
            neights[idx(i,j)*NEIGHCOUNT + 3] = (i==L-1) ? idx(0,j) : idx(i+1,j); //down
        }
    }
}

double eCalc(){
    int i;
    unsigned si, sj;
    eOld=0;
    mOld=0;
    mxOld=myOld=mzOld=0.;
    for (i=0;i<N;++i){
            si = parts[i];
            sj = parts[neights[i*NEIGHCOUNT + 0]];
            eOld += hamiltonianHeisenbergDMI(
                        pts[si*3+0], pts[si*3+1], pts[si*3+2],
                    pts[sj*3+0], pts[sj*3+1], pts[sj*3+2],
                    1, 0
                    );
            sj = parts[neights[i*NEIGHCOUNT + 1]];
            eOld += hamiltonianHeisenbergDMI(
                        pts[si*3+0], pts[si*3+1], pts[si*3+2],
                    pts[sj*3+0], pts[sj*3+1], pts[sj*3+2],
                    0, 1
                    );
            mxOld+=pts[si*3+0];
            myOld+=pts[si*3+1];
            mzOld+=pts[si*3+2];
    }
    mOld = sqrt(mxOld*mxOld+myOld*myOld+mzOld*mzOld);
    return eOld;
}

int main(int argc, char *argv[])
{
    (void)argc;(void)argv;

    rseed = 1;

    g = gsl_histogram2d_alloc(DOSSIZE1,DOSSIZE2);
    gsl_histogram2d_set_ranges_uniform(g,EMIN,EMAX,MMIN,MMAX);

    drop();
    genPoints();
    dumpPts();

    std::default_random_engine generator;
    generator.seed(rseed);
    std::uniform_int_distribution<ulong> distr(0,N-1);
    std::uniform_int_distribution<uint> uniform01(0, points);
    std::uniform_real_distribution<double> uniform02(0.0, 1.0);

    bool isFirstCycle=true, allowDump=false;
    unsigned accepted=0, rejected=0, total=0, cycles=0;
    double dx,dy,dz;
    unsigned newPoint, i, ss;
    size_t oldEIdx, newEIdx, oldMIdx, newMIdx;
    k=0;

    eCalc();
    gsl_histogram2d_find(g,eOld,mOld,&oldEIdx,&oldMIdx);


    for (i=0; i<DOSSIZE1*DOSSIZE2; ++i){
        g->bin[i]=0;
        h[i]=0;
        zeros[i]=0;
    }

    //g->bin[oldIdx]+=f;
    //h[oldIdx]+=1;
    //zeros[oldIdx]=true;

    //++total;

    cout<<"# einit = "<< eOld<<endl;
    dumpCfg(eOld);
    while (k<walks){
        //повторяем алгоритм сколько-то шагов
        for (i=0;i<N;i++){
            ulong partNum = distr(generator);

            newPoint = uniform01(generator);
            dx=pts[newPoint*3+0] - pts[parts[partNum]*3+0];
            dy=pts[newPoint*3+1] - pts[parts[partNum]*3+1];
            dz=pts[newPoint*3+2] - pts[parts[partNum]*3+2];

            de=0; dm=0;
            ss=parts[neights[partNum*NEIGHCOUNT + 0]];
            de += hamiltonianHeisenbergDMI(
                        dx, dy, dz,
                        pts[ss*3+0], pts[ss*3+1], pts[ss*3+2],
                        1, 0
                        );
            ss=parts[neights[partNum*NEIGHCOUNT + 1]];
            de += hamiltonianHeisenbergDMI(
                        dx, dy, dz,
                        pts[ss*3+0], pts[ss*3+1], pts[ss*3+2],
                        0, 1
                        );

            ss=parts[neights[partNum*NEIGHCOUNT + 2]];
            de += hamiltonianHeisenbergDMI(
                        pts[ss*3+0], pts[ss*3+1], pts[ss*3+2],
                        dx, dy, dz,
                        1, 0
                        );

            ss=parts[neights[partNum*NEIGHCOUNT + 3]];
            de += hamiltonianHeisenbergDMI(
                        pts[ss*3+0], pts[ss*3+1], pts[ss*3+2],
                        dx, dy, dz,
                        0, 1
                        );
            eNew=eOld+de;
            mNew = sqrt((mxOld+dx)*(mxOld+dx)+(myOld+dy)*(myOld+dy)+(mzOld+dz)*(mzOld+dz));

            gsl_histogram2d_find(g,eNew,mNew,&newEIdx,&newMIdx);

            if (uniform02(generator) <= exp(g->bin[oldEIdx*DOSSIZE2+oldMIdx]-g->bin[newEIdx*DOSSIZE2+newMIdx])) {
                eOld = eNew;
                mOld = mNew;
                mxOld += dx;
                myOld += dy;
                mzOld += dz;
                oldEIdx = newEIdx;
                oldMIdx = newMIdx;
                parts[partNum] = newPoint;
                ++accepted;
            } else {
                ++rejected;
            }

            g->bin[oldEIdx*DOSSIZE2+oldMIdx]+=f;
            h[oldEIdx*DOSSIZE2+oldMIdx]+=1;
            zeros[oldEIdx*DOSSIZE2+oldMIdx]=true;
            ++total;

            if (total % 1000000 == 0){ // renew energy every 100000 step for double error correction
                //double mxbkp=mxOld,mybkp=myOld,mzbkp=mzOld,mbkp=mOld,ebkp=eOld;
                eCalc();
                /*cout<<setprecision(5)<<total<<"\t"<<accepted<<"\t"<<rejected<<endl;
                cout<<mxOld<<"\t"<<myOld<<"\t"<<mzOld<<"\t"<<mOld<<"\t"<<eOld<<endl;
                cout<<mxbkp<<"\t"<<mybkp<<"\t"<<mzbkp<<"\t"<<mbkp<<"\t"<<ebkp<<endl;
                cout<<setprecision(17)<<(mxOld-mxbkp)<<"\t"<<myOld-mybkp<<"\t"<<mzOld-mzbkp<<"\t"<<mOld-mbkp<<"\t"<<eOld-ebkp<<endl<<endl;
                mxOld=mxbkp; myOld=mybkp; mzOld=mzbkp; mOld=mbkp; eOld=ebkp;*/
                gsl_histogram2d_find(g,eOld,mOld,&oldEIdx,&oldMIdx);
            }
        }
        ++cycles;

        //if (cycles%100==0){
            //cout<<"# total: "<<total<<"; accepted: "<<accepted<<"; rejected: "<<rejected<<"; f: "<<f<<endl;
        //}

        if (cycles>100000)
            isFirstCycle=false;

        //проверяем ровность диаграммы
        if (!isFirstCycle && cycles%10000==0){
            if (allowDump){
                dumpCfg(eOld);
                dump();
                cout<<"# dumped"<<endl;
                allowDump=false;
            }
            if (isFlat()){
                f/=2.;k++;
                for (i=0; i<DOSSIZE1*DOSSIZE2; ++i)
                    h[i]=0;
                norm();
                dump();
                //cout<<"accepted "<<(int)accepted<<endl;
                //cout<<"rejected "<<(int)rejected<<endl;
                cout<<"# step: "<<k<<"; h is flat, new f is "<<f<<endl;
                accepted=0; rejected=0; total=0; isFirstCycle = true;
            }
        }
    }


    dump();
}
