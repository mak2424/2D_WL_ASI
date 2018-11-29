#include <iostream>
#include <cstdlib>
#include <cstdbool>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <map>
#include <random>

#include <stdio.h>

#include <gsl/gsl_histogram2d.h>

using namespace std;

int PRECISION = 1e4;            // Точность 1eX, где X - Сколько знаков учитывать в энергии после запятой
                                // (1e0 - 0 знаков после запятой (для модели Изинга), 1e100 - 100 знаков после запятой)

unsigned N;              // number of spins in the lattice

signed char *spins;             //массив направления спинов. По умолчанию +1. n - число считанных спинов (число незакомментированных строк в csv-файле).
unsigned short *a_neighbours;   //число соседей каждого спина. Считается как число энергий в соответствующей строке в csv-файле.
unsigned short *neighbours;     //соседи каждого спина
unsigned int *sequencies;       //для каждого спина описывает, с какого ключа в массиве energies[] начинают описываться парные энергии
double *energies;               //сами энергии из файла. Описывается как одномерный массив. Длина массива - число парных энергий в csv-файле.
double *mx;                     //массив магнитных моментов спинов mx
double *my;                     //массив магнитных моментов спинов my
double emin, emax;              //минимумы и максимумы энергии
unsigned eCount=0;              //число пар энергий
unsigned long DOSSIZE1;         // number of bins in dos histogramm

int mmax, mmin;                 //минимумы и максимумы намагниченности
unsigned long DOSSIZE2;         // number of bins in dos histogramm

double accuracy=0.8;
double f=1;
unsigned walks=22;
int rseed;
unsigned k;

gsl_histogram2d *g;
unsigned *h;
char *zeros;

double  eNew, eOld,
        mNew, mOld,
        mxOld, mxNew,
        myOld, myNew;


/// Функция чтения файла с энергиями
int readCSV(char *filename){

    char c;                         //считанный из файла символ
    char symb[100000];                 //символ энергии в текстовом файле

    //get system sizes
    bool isFirstLine=true;
    N=0;
    FILE *file = fopen(filename, "r");

    if (!file)
        return 0;

    int fpos = 1, lastFpos=0;
    unsigned count_n=0;

    for(c=fgetc(file);c=='\n'||c=='\r'||c=='#';c=fgetc(file)){ //пропуск комментариев
        fgets(symb,100000,file);
    }

    fseek(file,-1,SEEK_CUR);       // сдвиг курсора на один символ назад
    int coursor=ftell(file);       // положение курсора начала данных

    do{
        c = fgetc(file);
        //        while(c=='#'){
        //            do c = fgetc(file2); while (c != '\n');           // нет необходимости, только если у нас не будет комментариев прямо посреди данных, но можно оставить
        //            c = fgetc(file2);
        //        }
        if (isFirstLine && c==';')
            ++N;
        if (c=='\n'){
            isFirstLine=false;
            count_n++;
        }

        if (c==';' || c=='\n') {
            if (fpos-1 != lastFpos)
                ++eCount;
            lastFpos = fpos;
        }

        fpos++;
    } while (c != EOF);
    ++N;
    if(count_n!=N)
        printf("!!!ERROR with number of element: Number of elements in first line does not correspond with number of lines");

    // reserve memory for arrays
    spins=(signed char *) malloc(N*sizeof(signed char));
    a_neighbours=(unsigned short *) malloc(N*sizeof(unsigned short));
    neighbours=(unsigned short *) malloc(eCount*sizeof(unsigned short));
    sequencies=(unsigned int *) malloc(N*sizeof(unsigned int));
    energies = (double *) malloc(eCount*sizeof(double));


    // read data

    fseek(file,coursor,SEEK_SET);      //устанавливаем курсор в начало данных

    double parsedNumber;
    int numInSymb=0;
    symb[0]='\0';
    int row=0;                  //line number in file (not account the commented lines)
    int col=0;                  //column number in line (taking to accound the ';' symbols)
    int neighCount=0;           //
    int energyNum=0;            //holds actual count of previously parsed energies
    eOld = 0;
    emax = 0;                   // сумма всех взаимодействий с положительным знаком
    long long eTemp=0;

    do {
        c = fgetc(file);

        if (c==';' || c=='\n' || c == EOF){ //if we found a number, process it
            if (numInSymb!=0){
                sscanf( symb, "%lf", &parsedNumber );
                neighbours[energyNum] = col;
                
                energies[energyNum] = parsedNumber;
                eTemp+=round(parsedNumber*10000000);
                emax += fabs(parsedNumber);

                numInSymb=0;
                ++neighCount;
                ++energyNum;
            }
            ++col;
        } else {
            if(c=='\r'){}
            else{
            symb[numInSymb] = c;
            symb[numInSymb+1] = '\0';
            ++numInSymb;
            }
        }

        if (c=='\n'){
            a_neighbours[row] = neighCount;
            sequencies[row] = energyNum-neighCount;
            col=0;
            neighCount=0;
            spins[row]=1;
            ++row;
        }
    } while (c != EOF);

    emax/=2;
    eOld=eTemp/20000000.;
    emax= ((double)ceil(emax * PRECISION)) / PRECISION; // округляем значение максимальной энергии до знака точности, в большую сторону
    emin = -emax;

    fclose(file);

    return 1;
}



/// Функция чтения файла с магнитными моментами
int read_mx_my(char *filename){

    char c;                         //считанный из файла символ
    char symb[100000];                 //символ магнитного момента в текстовом файле

    //get array sizes
    FILE *file = fopen(filename, "r");

    if (!file)
        return 0;

    unsigned count_n=0; //количество строк (элементов)

    for(c=fgetc(file);c=='\n'||c=='\r'||c=='#';c=fgetc(file)){ //пропуск комментариев
        fgets(symb,100000,file);
    }

    fseek(file,-1,SEEK_CUR);       // сдвиг курсора на один символ назад
    int coursor=ftell(file);       // положение курсора начала данных

    do{
        c = fgetc(file);
        //        while(c=='#'){
        //            do c = fgetc(file2); while (c != '\n');           // нет необходимости, только если у нас не будет комментариев прямо посреди данных, но можно оставить
        //            c = fgetc(file2);
        //        }
        if (c=='\n'){
            count_n++;
        }
    } while (c != EOF);
    
    // reserve memory for arrays
    mx = (double *) malloc(count_n*sizeof(double));
    my = (double *) malloc(count_n*sizeof(double));


    // read data

    fseek(file,coursor,SEEK_SET);      //устанавливаем курсор в начало данных

    double parsedNumber;
    int numInSymb=0;
    symb[0]='\0';
    int row=0;                  //line number in file (not account the commented lines)
    int col=0;                  //column number in line (taking to accound the ';' symbols)

    do {
        c = fgetc(file);

        if (c==';' || c=='\n' || c == EOF){ //if we found a number, process it
            if (numInSymb!=0){
                sscanf( symb, "%lf", &parsedNumber );
                
                if(col == 0)
                    mx[row] = parsedNumber;
                else if (col ==1)
                    my[row] = parsedNumber;

                numInSymb=0;
            }
            ++col;
        } else {
            if(c=='\r'){}
            else{
            symb[numInSymb] = c;
            symb[numInSymb+1] = '\0';
            ++numInSymb;
            }
        }

        if (c=='\n'){
            col=0;
            ++row;
        }
    } while (c != EOF);

    fclose(file);

    return 1;
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
    ofstream file(fname.str());
    file<<"# WL calculations for 2D spins on flat square lattice"<<endl;
    file<<"# with dipolar interaction"<<endl;
    file<<"# random seed="<<rseed<<endl;
    file<<"# bins number="<<DOSSIZE1<<" , "<<DOSSIZE2<<endl;
    file<<"# emin , emax = "<<emin<<" , "<<emax<<endl;
    file<<"# mmin , mmax = "<<mmin<<" , "<<mmax<<endl;

    for (unsigned long num1=0; num1<DOSSIZE1; ++num1){
        for (unsigned long num2=0; num2<DOSSIZE2; ++num2){
            if (zeros[num1*DOSSIZE2+num2]){
                file<<
                    num1<<"\t"<<
                    num2<<"\t"<<
                    g->xrange[num1]<<"\t"<<
                    g->yrange[num2]<<"\t"<<
                    g->bin[num1*DOSSIZE2+num2]<<"\t"<<
                    h[num1*DOSSIZE2+num2]<<endl;
            }
        }
    }
    file.close();
}

double recalcE(){
    unsigned i,j,is,js;
    eOld=0;
    mOld=0;
    mxOld=myOld=0.;
    for (i=0; i<N; ++i){
        for (j=sequencies[i]; j<sequencies[i]+a_neighbours[i]; ++j){
            is=spins[i];
            js=spins[neighbours[j]];
            if (is!=js)
                eOld-=energies[j];
            else
                eOld+=energies[j];
        }
        mxOld+=mx[i];
        myOld+=my[i];
    }
    eOld/=2.;
    mOld = sqrt(mxOld*mxOld+myOld*myOld);
    return eOld;
}


/// Переворот спина, подсчет изменения энергии
void rotate(int spin)
{
    unsigned i;
    long long dE=round(eOld*10000000);
    spins[spin] *= -1;
    for(i = sequencies[spin]; i<sequencies[spin]+a_neighbours[spin]; ++i){
        dE += round(energies[i]*spins[neighbours[i]]*spins[spin]*2*10000000);
    }
    eNew = dE/10000000.;

    mx[spin] *= -1;
    my[spin] *= -1;

    mxNew = mxOld+2*mx[spin];
    myNew = myOld+2*my[spin];
    mNew = sqrt(mxNew*mxNew + myNew*myNew);
}

/// Очистка памяти
void complete()
{
    // clean arrays
    free(spins);
    free(a_neighbours);
    free(neighbours);
    free(sequencies);
    free(energies);

    free(mx);
    free(my);

    free(g);
    free(h);
    free(zeros);
}


int main(int argc, char *argv[])
{
    (void)argc;(void)argv;

    rseed = 1;

    readCSV("new_csv.csv");
    read_mx_my("mx_my.csv");

    DOSSIZE1 = (int)((emax-emin)*PRECISION)+1;

    mmax = (N+1);
    mmin = (-N-1);
    DOSSIZE2 = mmax-mmin;

    h = (unsigned *) malloc(DOSSIZE1*DOSSIZE2*sizeof(unsigned));
    zeros = (char *) malloc(DOSSIZE1*DOSSIZE2*sizeof(char));
 
    cout << "emin = " << emin << "\nemax = " << emax << endl;
    cout << "mmin = " << mmin << "\nmmax = " << mmax << endl;
    cout << "DOSSIZE1 = " << DOSSIZE1 << "\nDOSSIZE2 = " << DOSSIZE2 << endl;

    g = gsl_histogram2d_alloc(DOSSIZE1,DOSSIZE2);
    gsl_histogram2d_set_ranges_uniform(g,emin,emax,mmin,mmax);

    std::default_random_engine generator;
    generator.seed(rseed);
    std::uniform_int_distribution<unsigned long> distr(0,N-1);
    std::uniform_real_distribution<double> uniform02(0.0, 1.0);

    bool isFirstCycle=true, allowDump=false;
    unsigned accepted=0, rejected=0, total=0, cycles=0;
    size_t oldEIdx, newEIdx, oldMIdx, newMIdx;
    k=0;

    recalcE();
    gsl_histogram2d_find(g,eOld,mOld,&oldEIdx,&oldMIdx);

    for (unsigned i=0; i<DOSSIZE1*DOSSIZE2; ++i){
        g->bin[i]=0;
        h[i]=0;
        zeros[i]=0;
    }

    cout<<"# einit = "<< eOld<<endl;

    while (k<walks){
        //повторяем алгоритм сколько-то шагов
        for (unsigned i=0;i<N;i++){
            unsigned long partNum = distr(generator);//выбираем случайный спин

            rotate(partNum);//получаем eNew, mNew

            gsl_histogram2d_find(g,eNew,mNew,&newEIdx,&newMIdx);

            if (uniform02(generator) <= exp(g->bin[oldEIdx*DOSSIZE2+oldMIdx]-g->bin[newEIdx*DOSSIZE2+newMIdx])) {
                eOld = eNew;
                mOld = mNew;
                mxOld = mxNew;
                myOld = myNew;

                oldEIdx = newEIdx;
                oldMIdx = newMIdx;
                ++accepted;
            } else {
                spins[partNum] *= -1;
                ++rejected;
            }

            g->bin[oldEIdx*DOSSIZE2+oldMIdx]+=f;
            h[oldEIdx*DOSSIZE2+oldMIdx]+=1;
            zeros[oldEIdx*DOSSIZE2+oldMIdx]=true;
            ++total;

            if (total % 1000000 == 0){ // renew energy every 100000 step for double error correction
                //double mxbkp=mxOld,mybkp=myOld,mbkp=mOld,ebkp=eOld;
                recalcE();
                /*cout<<setprecision(5)<<total<<"\t"<<accepted<<"\t"<<rejected<<endl;
                cout<<mxOld<<"\t"<<myOld<<"\t"<<"\t"<<mOld<<"\t"<<eOld<<endl;
                cout<<mxbkp<<"\t"<<mybkp<<"\t"<<"\t"<<mbkp<<"\t"<<ebkp<<endl;
                cout<<setprecision(17)<<(mxOld-mxbkp)<<"\t"<<myOld-mybkp<<"\t"<<"\t"<<mOld-mbkp<<"\t"<<eOld-ebkp<<endl<<endl;
                mxOld=mxbkp; myOld=mybkp; mOld=mbkp; eOld=ebkp;*/
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
                dump();
                cout<<"# dumped"<<endl;
                allowDump=false;
            }
            if (isFlat()){
                f/=2.;k++;
                for (unsigned i=0; i<DOSSIZE1*DOSSIZE2; ++i)
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
    
    complete(); //очистка памяти
}
