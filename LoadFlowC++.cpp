//TO PREPARE THE PROGRAM, THE GLOBAL VARIABLE Y's ARRAY DIMENSIONS
//MUST BE FITTED TO THE APPROPRIATE VALUES (IT'S A SQUARE MATRIX
//WITH DIMENSION EQUAL TO THE # OF SYSTEM BUSES). THE YBUS.TXT FILE
//CONTAINS AN INITIAL ENTRY CONTAINING THE # OF BUSES N, FOLLOWED BY A
//SEQUENTIAL LIST OF REAL PARTS OF THE YBUS ENTRIES. THE ABUTTING LIST
//AFTERWARDS (FROM LIST ELEMENT N^2+1) CONTAINS THE IMAGINARY YBUS ENTRIES
//N is set to the number of system buses

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <complex>
#include <iostream>

using namespace std;
typedef complex<double> com; //declare complex variables with "com"

double h=0.000000001;
int MaxNoIterations=25;
double MinimumNormV=0.0000001;
int N=14 ; //# of buses
com Y[15][15]; //dimension of Ybus +1 (we want to reference our array
               //by its corresponding matrix row/column indices; but array
double sysbase=100;     //MVA basis of the system



double* inv(double**J, int nv, double dP[]);
void LoadYbus();
void LoadConditions();
bool NewtonRaphson(com Vfull[],double V[],double Ss[],int nv,char Btype[],double Vhold[]);
com* Powers(com Vfull[]);



int main(){

    LoadYbus();     // read all the Ybus entries from the text file Ybus.txt

    FILE * pFile; //pFile contains the address of a file object
   pFile = fopen ("OperatingConditions.txt","r"); //Everything from here until the Newton Raphson iterations is to load the operating conditions of the system (scheduled powers, initial guesses)

       int c; char Btype[N+1]; double Vdata[N+1][3]; //variable 1 - voltage magnitude
                                                 //            2 - delta
               double Sdata [N+1][3]; //secondary index  0 - real generation on bus n
                              //                        1 - real part of Pscheduled
                              //                        2 - imag part of Pscheduled

                int n; int npv=0; double temp; char tempc; int tempi;
    while (! feof (pFile)) {
     fscanf(pFile, "%d", &tempi); n=tempi;                  //bus number
       fscanf(pFile, "%c", &tempc); Btype [n] = tempc;      //type of bus
          if (tempc=='v') {npv++;}    //count the number of PV buses (needed to determine the number of variables in the system)
       fscanf(pFile, "%lf", &temp); Vdata[n][1]=temp;        //voltage magnitude
       fscanf(pFile, "%lf", &temp); Vdata[n][2]=temp;        //voltage angle
       fscanf(pFile, "%lf", &temp); Sdata[n][0]=temp/sysbase;    //generated/injected real power
       fscanf(pFile, "%lf", &temp); Sdata[n][1]=temp/sysbase;    //dissipated/demanded real power
       fscanf(pFile, "%lf", &temp); Sdata[n][2]=temp/sysbase;    //dissipated/demanded complex power (for complex power generation, subtract it from this value in the OperatingConditions text file)

                                                }
     fclose (pFile);

int nv=2*(N-1)-npv;     //number of variables (voltages + angles) = 2*(#buses-1)-(#PV buses)
    char Vtype[nv+1]; char Stype[nv+1];  //stores the type of each variable (voltage magnitude/angle, real power/complex power)
    double V[nv+1];
    double Ss[nv+1]; int varcount=1; //create the variables column matrix, Sscheduled, and index the current variable in question with varcount

 for (int b=1;b<=n;b++){

  if (Btype[b] == 's') {} //if bus b is the slack bus

    if (Btype[b] == 'q') {   //if bus b is a PQ bus (load the voltage 'v' and angle 'd' into the variables column matrix V
        V[varcount]=Vdata[b][1]; V[varcount+1]=Vdata[b][2];
          Vtype [varcount]='v'; Vtype[varcount+1]='d';  //"v" denotes magnitude, "d" angle
        Ss[varcount]=Sdata[b][1]-Sdata[b][0]; Ss[varcount+1]=Sdata[b][2];
          Stype[varcount]='p'; Stype[varcount+1]='q'; //"p" denotes real power, "q" imag. power
            varcount+=2;
    }
    if (Btype[b] == 'v') {  //if bus b is a PV bus (load only the angle 'd' into the variables column matrix V
        V[varcount]=Vdata[b][2];
           Vtype [varcount]='d';    //this variable is of type delta (angle)
        Ss[varcount]=Sdata[b][1]-Sdata[b][0];   //scheduled power on a PV bus is Sdemanded-Sgenerated
           Stype[varcount]='q';     //this power is of type q (complex power)
              varcount++;
    }
}

    com Vfull[N+1];
   for (int a=1;a<=N;a++){
      Vfull[a]=polar(Vdata[a][1],Vdata[a][2]);
    }   //Vfull contains voltage magnitude and phase information for ALL buses (used to find unperturbed and perturbed powers)

        double* Vhold=(double *)malloc((nv+1)*sizeof(double));
        for(int i=1;i<=nv;i++) Vhold[i]=0;

          //NEWTON RAPHSON ITERATIONS BEGIN HERE
            bool keepgoing=1;
            int iteration=1;
        while (keepgoing and iteration<MaxNoIterations+1) {iteration++;
        keepgoing=NewtonRaphson(Vfull,V,Ss,nv,Btype,Vhold);} //perturb variables, find jacobian, use LU decomposition to find the dV (change in voltages) matrix, return the new "guess" (=old guess + dV)

    //system("pause");
   return 0;
}


bool NewtonRaphson(com Vfull[],double V[],double Ss[],int nv,char Btype[],double Vhold[]){

  static int iteration=1;    //this is incremented at the bottom of this function, and retains its old value in successive calls to this function
  if (iteration>1) {for(int i=1;i<=nv;i++) {V[i]=Vhold[i];}}

         int n=1;
         int v=1;
    while (v<=nv and n<=N){   //create a column matrix Vfull with the voltage and angle of EVERY bus using the variables column matrix
            if      (Btype[n]=='q') {Vfull[n]=polar(V[v],V[v+1]); n++; v+=2;}
            else if (Btype[n]=='v') {Vfull[n]=polar(abs(Vfull[n]),V[v]); n++; v++;}
            else if (Btype[n]=='s') {n++;}
    }

 //printf("(%6.4f, %6.4f)\n", Vfull[c]);

    com* Sgfull=Powers(Vfull); //starts from 0 (these are the powers for each bus at our voltage guess)
    double Sg[nv+1]; //starts from 1

             n=0;
             v=1;
        while (v<=nv and n<N){   //Put powers at our guess Sgfull into the Sg variables column matrix for the Jacobian
            if      (Btype[n+1]=='q') {Sg[v]=Sgfull[n].real(); Sg[v+1]=Sgfull[n].imag(); n++; v+=2;}
            else if (Btype[n+1]=='v') {Sg[v]=Sgfull[n].real(); n++; v++;}
            else if (Btype[n+1]=='s') {n++;}
        }

        double dP[nv+1];
         for (int i=1;i<=nv;i++) dP[i]=Sg[i]-Ss[i];   //delta powers




      double **J=(double **)malloc((nv+1)*sizeof(double *));      //allocate memory for the two-dimensional Jacobian, so that it can be accessed outside our current function in the LU decomposition implementing function (inv)
      for(int i=1;i<=nv;i++)J[i]=(double *)malloc((nv+1)*sizeof(double));
                for(int r=1;r<=nv;r++)for(int c=1;c<=nv;c++) J[r][c]=0;

                com Vper[N+1];  //perturb variables, find S, set up each column of the Jacobian
                for (int v=1;v<=N;v++) Vper[v]=Vfull[v]; //initialize perturbed variables array to the full set of unperturbed voltages and angles
                int vp=1;  //indexes the variable being perturbed
                com* Sp;

            for (int n=1;n<=N;n++){  //indexes the bus from which the variable being perturbed is derived from the bus type (Btype)
                    if      (Btype[n]=='s') {}  //if the bus in question is the slack bus, there is nothing to load into the jacobian from it
                    else if (Btype[n]=='q') {     // if the bus in question is a PQ bus, we need to perturb a voltage and an angle, find the corresponding powers, and fill two jacobian columns with this information

                        Vper[n]=polar(abs(Vfull[n])+h,arg(Vfull[n])); //first perturb the voltage magnitude at this bus
                        Sp=Powers(Vper); //find powers with this perturbed magnitude

                                    int jr=1; //jacobian row (jacobian column=vp)
                                for (int jb=0;jb<N;jb++){ //index to the bus for the jacobian
                                    if      (Btype[jb+1]=='s'){}
                                    else if (Btype[jb+1]=='q'){J[jr][vp]=(Sgfull[jb].real()-Sp[jb].real())/h; jr++;  //use h even though voltage/delta - P/Q are in different coordinates?
                                                               J[jr][vp]=(Sgfull[jb].imag()-Sp[jb].imag())/h; jr++;}
                                    else if (Btype[jb+1]=='v'){J[jr][vp]=(Sgfull[jb].real()-Sp[jb].real())/h; jr++;}
                                }


                                vp++; // move onto perturbing the next variable on the PQ bus
                                   Vper[n]=polar(abs(Vfull[n]),arg(Vfull[n])+h); //reset our last perturbation and perturb the angle now
                                Sp=Powers(Vper); //find powers with this perturbed angle

                                    jr=1; //jacobian row (jacobian column=vp)
                                for (int jb=0;jb<N;jb++){ //index to the bus for the jacobian
                                    if      (Btype[jb+1]=='s'){}
                                    else if (Btype[jb+1]=='q'){
                                            J[jr][vp]=(Sgfull[jb].real()-Sp[jb].real())/h; jr++;  //use h even though voltage/delta - P/Q are in different coordinates?
                                            J[jr][vp]=(Sgfull[jb].imag()-Sp[jb].imag())/h; jr++;
                                    }
                                    else if (Btype[jb+1]=='v'){J[jr][vp]=(Sgfull[jb].real()-Sp[jb].real())/h; jr++;}
                                }

                                Vper[n]=Vfull[n]; // reset the perturbation column
                                    vp++;   //proceed to perturb variables on the next bus
                    }

                    else if (Btype[n]=='v'){    //if the bus in question (for loading the jacobian) is a PV bus

                            Vper[n]=polar(abs(Vfull[n]),arg(Vfull[n])+h); //set up the perturbation column with a perturbed delta at bus n
                            Sp=Powers(Vper); //powers perturbed

                                int jr=1; //jacobian row (jacobian column = vp = the index to the perturbed variable)
                                for (int jb=0;jb<N;jb++){ //index to the bus for the jacobian
                                    if      (Btype[jb+1]=='s'){}
                                    else if (Btype[jb+1]=='q'){J[jr][vp]=(Sgfull[jb].real()-Sp[jb].real())/h; jr++;  //use h even though voltage/delta - P/Q are in different coordinates?
                                                               J[jr][vp]=(Sgfull[jb].imag()-Sp[jb].imag())/h; jr++;}
                                    else if (Btype[jb+1]=='v'){J[jr][vp]=(Sgfull[jb].real()-Sp[jb].real())/h; jr++;}
                                }
                                                        /*
                                                        static int abcd=0; if(abcd==0){
                                                        cout<<endl<<"V perturbed"<<endl<<endl;
                                                        for (int r=1;r<=N;r++) cout<<Vper[r]<<endl;} abcd++;
                                                        */
                                            Vper[n]=Vfull[n];
                                            vp++;   // proceed to perturb variables on the next bus
                    }
            }
/*
static int abcd=0; if(abcd==0){
cout<<endl<<"Jacobian column 1"<<endl<<endl;
     for (int r=1;r<=nv;r++) cout<<J[r][vp]<<endl; abcd++;
     cout<<endl<<"Sgfull"<<endl;
     for (int r=0;r<N;r++) cout<<Sgfull[r]<<endl;
     cout<<endl<<"Sp"<<endl;
     for (int r=0;r<N;r++) cout<<Sp[r]<<endl; */

//                cout<<endl<<"Jacobian"<<endl<<endl;
//      for (int r=1;r<=11;r++) {cout<<endl<<"column "<<r<<endl; for (int c=1;c<=nv;c++) cout<<J[c][r]<<endl;}

     double* dV=inv(J,nv,dP);

        //THE REST OF THE FUNCTION JUST PRINTS ALL VOLTAGES AND ANGLES WITH THE dV ADDED
        //CALCULATES NORMS, RETURN A BOOLEAN FOR THE CONVERGENCE CRITERIA DECISION

//    cout<<"dV matrix"<<endl;
//    for (int r=1;r<=nv;r++){
//    cout<<dV[r]<<endl; }

    for (int i=1;i<=nv;i++) V[i]+=dV[i];
            for (int i=1;i<=nv;i++) {Vhold[i]=V[i];} //cout<<endl<<"Vi "<<V[i]<<endl;} // put the voltage guess away for the next iteration
      n=1;
      v=1;
    while (v<=nv and n<=N){   //create a vector with real and complex voltage information for every bus using the updated variable values
            if      (Btype[n]=='q') {Vfull[n]=polar(V[v],V[v+1]); n++; v+=2;}
            else if (Btype[n]=='v') {Vfull[n]=polar(abs(Vfull[n]),V[v]); n++; v++;}
            else if (Btype[n]=='s') {n++;}
    }

    double normV=0; for (int i=1;i<=nv;i++) normV+=dV[i]*dV[i];
              normV=pow(normV,0.5);
    com* Sfinal=Powers(Vfull);

    cout<<endl<<"dP                         dV               Sguess               Ssch"<<endl;
         for (int i=1;i<=nv;i++) cout<<dP[i]<<"         "<<dV[i]<<"         "<<Sg[i]<<"         "<<Ss[i]<<endl;

    cout<<endl<<endl<<endl<<"Iteration "<<iteration<<endl;
                    iteration++;
    cout<<endl<<"Voltages/Powers (Bus#  Vmag  Vdelta  <Pnet , Qnet>)"<<endl<<endl;
    for (int i=1;i<=N;i++) cout<<i<<"        "<<abs(Vfull[i])<<"        "<<arg(Vfull[i])<<"        "<<Sfinal[i-1]<<endl;
    cout<<endl<<"Voltages norm = "<<normV<<endl;

       bool keepgoing=1;
      if (normV<MinimumNormV) {keepgoing=0;}
      return keepgoing;

}



com* Powers(com* Vfull){

                //com S[N+1];
        com* S=(com *)malloc(N*sizeof(com));    //allocate the memory needed to access the powers outside of the function
            for (int i=0;i<N;i++){S[i]=com(0,0);}


        com YV[N+1];
            for (int i=1;i<=N;i++){YV[i]=com(0,0);}

        //com* pS=S;
            for (int r=1;r<=N;r++){//YV[r]=com(0,0); // YV
            for (int c=1;c<=N;c++){YV[r]+=Y[r][c]*Vfull[c];}}

            for (int r=1;r<=N;r++) {YV[r]=com(YV[r].real(),-(YV[r].imag()));} //YV*

            for (int r=0;r<N;r++){//S[r]=com(0,0); // V(YV)* = S
                                   S[r]=Vfull[r+1]*YV[r+1];
                                   }
                return S;
}


void LoadConditions(){

    FILE * pFile; //pFile contains the address of a file object
   pFile = fopen ("OperatingConditions.txt","r");

       int c; char Btype[N+1]; double Vdata[N+1][3]; //variable 1 - voltage magnitude
                                                    // 2 - delta
               double Sdata [N+1][3]; // secondary index 0 - real generation on bus n
                              //                        1 - real part of Pscheduled
                              //                        2 - imag part of Pscheduled

                int n; int npv=0; double temp; char tempc; int tempi;
    while (! feof (pFile)) {
     fscanf(pFile, "%d", &tempi); n=tempi;
       fscanf(pFile, "%c", &tempc); Btype [n] = tempc;

          if (tempc=='v') {npv++;}
       fscanf(pFile, "%f", &temp); Vdata[n][1]=temp;
       fscanf(pFile, "%f", &temp); Vdata[n][2]=temp;
       fscanf(pFile, "%f", &temp); Sdata[n][0]=temp;
       fscanf(pFile, "%f", &temp); Sdata[n][1]=temp;
       fscanf(pFile, "%f", &temp); Sdata[n][2]=temp;

    }
     fclose (pFile);



int nv=2*(N-1)-npv;
    char Vtype[nv+1]; char Stype[nv+1];
    double V[nv+1];
    double S[nv+1]; int varcount=1;

 for (int b=1;b<=n;b++){

  if (Btype[b] == 's') {}

    if (Btype[b] == 'q') {
        V[varcount]=Vdata[b][1]; V[varcount+1]=Vdata[b][2];
          Vtype [varcount]='v'; Vtype[varcount+1]='d';  //"v" denotes magnitude, "d" angle
        S[varcount]=Sdata[b][1]-Sdata[b][0]; S[varcount+1]=Sdata[b][2];
          Stype[varcount]='p'; Stype[varcount+1]='q'; //"p" denotes real power, "q" imag. power
            varcount+=2;
    }
    if (Btype[b] == 'v') {
        V[varcount]=Vdata[b][2];
           Vtype [varcount]='d';
        S[varcount]=Sdata[b][1]-Sdata[b][0];
           Stype[varcount]='p';
              varcount++;
    }
 }

//            cout<<endl<<"V matrix inital guess"<<endl<<endl;
//                    for (int r=1;r<=nv;r++){
    //for (int c=1;c<=nv;c++){
//      cout<<V[r]<<endl; }//}
//      cout<<endl<<"S matrix"<<endl<<endl;
//      for (int r=1;r<=nv;r++){
    //for (int c=1;c<=nv;c++){
      //cout<<S[r]<<endl; }//}

//cout<<nv<<endl;



//         if ( fgets (buffer , 100 , pFile) != NULL )
//         fputs (buffer , stdout);
}



void LoadYbus(){
       FILE * pFile; //pFile contains the address of a file object
       pFile = fopen ("Ybus.txt","r");
        //fputs ("fopen example",pFile);
        //fclose (pFile);

       int yd;
       fscanf(pFile, "%d", &yd);
       if (yd == NULL) perror ("Error opening Ybus");
       else{   //cout<<yd<<" Bus System Ybus entries:"<<endl;
        //com Y[yd+1][yd+1];
        //com *pY=Y[1][1];
           int c=1; int r=1;
           double temp=0; //double rpart;
           //Y[1][0]=yd; Y[0][1]=yd;

         while ( ! feof (pFile) ){
           if ( fscanf(pFile, "%lf", &temp) != NULL ){
            if (r<=yd){ Y[r][c]=com(temp,1);}
            else {Y[r-yd][c]=com(Y[r-yd][c].real(),temp);}

             if   (c<yd) {c++;}
             else {c=1; r++;}
            }
    //         if ( fgets (buffer , 100 , pFile) != NULL )
    //         fputs (buffer , stdout);
         }
        }
     fclose (pFile);
     }


double* inv(double** J,int nv, double dP[]){

int rl,cl,ru,cu;
double backadd=0;
double Z[nv+1]; for (int c=1;c<=nv;c++) Z[c]=0;
double* dV=(double *)malloc((nv+1)*sizeof(double));   // allocate and initialize the change in variables matrix
    for (int i=1;i<=nv;i++) dV[i]=0;

double L[nv+1][nv+1]; double U[nv+1][nv+1];
for (int c=1;c<=nv;c++) for (int r=1;r<=nv;r++) L[r][c]=0;
for (int c=1;c<=nv;c++) for (int r=1;r<=nv;r++) U[r][c]=0;
for (int c=1;c<=nv;c++) L[c][c]=1;


    //DECOMPOSE INTO L AND U
for (cu=1;cu<=nv;cu++) for (rl=1;rl<=nv;rl++){

    if (rl<=cu) {
        for (cl=1;cl<rl;cl++){backadd+=L[rl][cl]*U[cl][cu];}
        U[rl][cu]=J[rl][cu]-backadd;
        //cout<<"entry "<<backadd<<endl;
    }
    else {
        for (cl=1;cl<rl;cl++){backadd+=L[rl][cl]*U[cl][cu];}
        L[rl][cu]=(J[rl][cu]-backadd)/U[cu][cu];
         if (U[cu][cu]==0) {L[rl][cu]=0;}
    }
         backadd=0;
}


    //FORWARD SUBSTITUTION IN L
for (int t=1;t<=nv;t++){
for (cl=1;cl<t;cl++){
backadd+=L[t][cl]*Z[cl];
                    }
Z[t]=dP[t]-backadd;
backadd=0;      }


    //BACK SUBSTITUTION IN U
for (int t=nv;t>0;t--){ for (cu=nv;cu>t;cu--){
backadd+=U[t][cu]*dV[cu];}

dV[t]=(Z[t]-backadd)/U[t][t];
backadd=0;
}

/*
                                    static int abcd=0; if(abcd==0){
cout<<endl<<"J matrix"<<endl<<endl;
     for (int r=1;r<=nv;r++)for (int c=1;c<=nv;c++) cout<<J[r][c]<<endl; abcd++;} */
//      cout<<endl<<"Sg full matrix"<<endl<<endl;
//      for (int r=1;r<=nv;r++) cout<<Sg[r]<<endl;

return dV;
}

