#include <iostream>
#include <cstdlib>
#include <memory>
#include <iomanip>
#include <cmath>
#include <complex>
#include "Complex.h"
#include <sys/time.h>

using namespace std;

#define nstart 0
#define nend 3

//Outputs are ssxa and scha, rest of the passed parameters are the inputs
inline void ssxt_scht_solver(double wxt, int igp, int my_igp, int ig, GPPComplex wtilde, GPPComplex wtilde2, GPPComplex Omega2, GPPComplex matngmatmgp, GPPComplex matngpmatmg, GPPComplex mygpvar1, GPPComplex mygpvar2, GPPComplex& ssxa, GPPComplex& scha, GPPComplex I_eps_array_igp_myIgp)
{
    GPPComplex expr0( 0.0 , 0.0);
    double ssxcutoff;
    double to1 = 1e-6;
    double sexcut = 4.0;
    double limitone = 1.0/(to1*4.0);
    double limittwo = 0.5*0.5;
    GPPComplex sch(0.00, 0.00);
    GPPComplex ssx(0.00, 0.00);

    GPPComplex wdiff = doubleMinusGPPComplex(wxt , wtilde);

    GPPComplex cden = wdiff;
    double rden = 1/GPPComplex_real(GPPComplex_product(cden , GPPComplex_conj(cden)));
    GPPComplex delw = GPPComplex_mult(GPPComplex_product(wtilde , GPPComplex_conj(cden)) , rden);
    double delwr = GPPComplex_real(GPPComplex_product(delw , GPPComplex_conj(delw)));
    double wdiffr = GPPComplex_real(GPPComplex_product(wdiff , GPPComplex_conj(wdiff)));

    if((wdiffr > limittwo) && (delwr < limitone))
    {
        sch = GPPComplex_product(delw , I_eps_array_igp_myIgp);
        double cden = wxt*wxt;
        rden = cden*cden;
        rden = 1.00 / rden;
        ssx = GPPComplex_mult(Omega2 , cden , rden);
    }
    else if (delwr > to1)
    {
        sch = expr0;
        cden = GPPComplex_mult(GPPComplex_product(wtilde2, doublePlusGPPComplex((double)0.50, delw)), 4.00);
//        cden = (double) 4.00 * wtilde2 * (delw + (double)0.50);
        rden = GPPComplex_real(GPPComplex_product(cden , GPPComplex_conj(cden)));
        rden = 1.00/rden;
        ssx = GPPComplex_product(GPPComplex_product(-Omega2 , GPPComplex_conj(cden)), GPPComplex_mult(delw, rden));
    }
    else
    {
        sch = expr0;
        ssx = expr0;
    }

    ssxcutoff = GPPComplex_abs(I_eps_array_igp_myIgp) * sexcut;
    if((GPPComplex_abs(ssx) > ssxcutoff) && (wxt < 0.00)) ssx = expr0;

    ssxa += GPPComplex_product(matngmatmgp , ssx);
    scha += GPPComplex_product(matngmatmgp , sch);
}

//This function writes its results to achstemp, rest of the parameters are its inputs.
inline void reduce_achstemp(int n1, int* inv_igp_index, int ncouls, GPPComplex *aqsmtemp, GPPComplex *aqsntemp, GPPComplex *I_eps_array, GPPComplex& achstemp, int ngpown, double* vcoul)
{
    double to1 = 1e-6;
    for(int my_igp = 0; my_igp< ngpown; my_igp++)
    {
        GPPComplex schstemp(0.0, 0.0);
        GPPComplex schs(0.0, 0.0);
        GPPComplex matngmatmgp(0.0, 0.0);
        GPPComplex matngpmatmg(0.0, 0.0);
        GPPComplex mygpvar1(0.00, 0.00);
        GPPComplex mygpvar2(0.00, 0.00);
        int igp = inv_igp_index[my_igp];
        if(igp >= ncouls)
            igp = ncouls-1;

        if(!(igp > ncouls || igp < 0))
        {

            mygpvar1 = GPPComplex_conj(aqsmtemp[n1*ncouls+igp]);
            mygpvar2 = aqsntemp[n1*ncouls+igp];
            GPPComplex schs_pos = I_eps_array[my_igp*ncouls+igp];

            schs = -schs_pos;
            matngmatmgp = GPPComplex_product(aqsntemp[n1*ncouls+igp] , mygpvar1);

            if(GPPComplex_abs(schs) > to1)
                GPPComplex_fma(schstemp, matngmatmgp, schs);
        }
        else
        {
            for(int ig=1; ig<ncouls; ++ig)
                GPPComplex_fms(schstemp, GPPComplex_product(aqsntemp[n1*ncouls+ig], I_eps_array[my_igp*ncouls+ig]), mygpvar1);
        }
        achstemp += GPPComplex_mult(schstemp , vcoul[igp] *(double) 0.5);
    }
}



//Performs the calculation for the first nvband iterations.
//Outputs are ssxt and scht, rest of the passed parameters are the inputs
inline void flagOCC_solver(double wxt, GPPComplex *wtilde_array, int my_igp, int n1, GPPComplex *aqsmtemp, GPPComplex *aqsntemp, GPPComplex *I_eps_array, GPPComplex &ssxt, GPPComplex &scht, int ncouls, int igp, GPPComplex *ssxa, GPPComplex* scha)
{
    GPPComplex matngmatmgp = GPPComplex(0.0, 0.0);
    GPPComplex matngpmatmg = GPPComplex(0.0, 0.0);
    GPPComplex expr0(0.00, 0.00);
    for(int ig=0; ig<ncouls; ++ig)
    {
        ssxa[ig] = expr0;
        scha[ig] = expr0;
        GPPComplex wtilde = wtilde_array[my_igp*ncouls+ig];
        GPPComplex wtilde2 = GPPComplex_square(wtilde);
        GPPComplex Omega2 = GPPComplex_product(wtilde2,I_eps_array[my_igp*ncouls+ig]);
        GPPComplex mygpvar1 = GPPComplex_conj(aqsmtemp[n1*ncouls+igp]);
        GPPComplex mygpvar2 = aqsmtemp[n1*ncouls+igp];
        GPPComplex matngmatmgp = GPPComplex_product(aqsntemp[n1*ncouls+ig] , mygpvar1);
        if(ig != igp) matngpmatmg = GPPComplex_product(GPPComplex_conj(aqsmtemp[n1*ncouls+ig]) , mygpvar2);

        ssxt_scht_solver(wxt, igp, my_igp, ig, wtilde, wtilde2, Omega2, matngmatmgp, matngpmatmg, mygpvar1, mygpvar2, ssxa[ig], scha[ig], I_eps_array[my_igp*ncouls+ig]); 
        ssxt += ssxa[ig];
        scht += scha[ig];
    }
}

//Outputs is scht, rest of the passed parameters are the inputs
inline void noflagOCC_solver(double wxt, GPPComplex *wtilde_array, int my_igp, int n1, GPPComplex *aqsmtemp, GPPComplex *aqsntemp, GPPComplex *I_eps_array, GPPComplex &ssxt, GPPComplex &scht, int ncouls, int igp, GPPComplex *scha)
{
    double limittwo = 0.5*0.5;
    GPPComplex mygpvar1 = GPPComplex_conj(aqsmtemp[n1*ncouls+igp]);
    GPPComplex scht_loc(0.00, 0.00);
    
    for(int ig = 0; ig<ncouls; ++ig)
    {
        GPPComplex wdiff = doubleMinusGPPComplex(wxt , wtilde_array[my_igp*ncouls+ig]);
        double wdiffr = GPPComplex_real(GPPComplex_product(wdiff , GPPComplex_conj(wdiff)));
        double rden = 1/wdiffr;

        GPPComplex delw = GPPComplex_mult(GPPComplex_product(wtilde_array[my_igp*ncouls+ig] , GPPComplex_conj(wdiff)) ,rden); 
        double delwr = GPPComplex_real(GPPComplex_product(delw , GPPComplex_conj(delw)));

        scht_loc += GPPComplex_product(GPPComplex_product(mygpvar1 , aqsntemp[n1*ncouls+ig]) , GPPComplex_product(delw , I_eps_array[my_igp*ncouls+ig])) ;
    }

    scht = scht_loc;
}

//This function calculates the first nvband iterations of the outermost loop
//output achtemp, asxtemp, acht_n1_loc
inline void till_nvband(int nvband, int ngpown, int ncouls, int *inv_igp_index, double *wx_array, GPPComplex *wtilde_array, GPPComplex *aqsmtemp, GPPComplex *aqsntemp, GPPComplex *I_eps_array, GPPComplex *asxtemp, double *vcoul, GPPComplex *acht_n1_loc, GPPComplex *achtemp, GPPComplex *ssx_array, GPPComplex *sch_array, const double occ, GPPComplex *ssxa, GPPComplex *scha)
{
    GPPComplex scht(0.00, 0.00);
    GPPComplex ssxt(0.00, 0.00);
    GPPComplex expr0(0.00, 0.00);
    for(int n1 = 0; n1<nvband; ++n1) 
    {
        for(int my_igp=0; my_igp<ngpown; ++my_igp)
        {
            int igp = inv_igp_index[my_igp];
            if(igp >= ncouls)
                igp = ncouls-1;

            for(int iw=nstart; iw<nend; iw++)
            {
                ssx_array[iw] = expr0;
                sch_array[iw] = expr0;
                scht = ssxt = expr0;
                double wxt = wx_array[iw];
                flagOCC_solver(wxt, wtilde_array, my_igp, n1, aqsmtemp, aqsntemp, I_eps_array, ssxt, scht, ncouls, igp, ssxa, scha);
                ssx_array[iw] += ssxt;
                sch_array[iw] +=GPPComplex_mult(scht, 0.5);
                
                asxtemp[iw] += GPPComplex_mult(ssx_array[iw] , occ * vcoul[igp]);//Store output of the first nvband iterations.
            }

            for(int iw=nstart; iw<nend; ++iw)
                achtemp[iw] += GPPComplex_mult(sch_array[iw] , vcoul[igp]);//Store final output here

            acht_n1_loc[n1] += GPPComplex_mult(sch_array[2] , vcoul[igp]);
        }
    }
}



inline void noFlagOCCSolver(int n1, int nvband, int ngpown, int ncouls, int *inv_igp_index, double *wx_array, GPPComplex *wtilde_array, GPPComplex *aqsmtemp, GPPComplex *aqsntemp, GPPComplex *I_eps_array, double *vcoul, GPPComplex *achtemp, GPPComplex *sch_array, GPPComplex *acht_n1_loc, GPPComplex *scha)
{
    GPPComplex expr0(0.00, 0.00);
    GPPComplex ssxt(0.00, 0.00);
    GPPComplex scht(0.00, 0.00);

    for(int my_igp=0; my_igp<ngpown; ++my_igp)
    {
        int igp = inv_igp_index[my_igp];
        if(igp >= ncouls)
            igp = ncouls-1;

        for(int iw=nstart; iw<nend; ++iw)
        {
            sch_array[iw] = expr0;
            scht = ssxt = expr0;
            double wxt = wx_array[iw];
            noflagOCC_solver(wxt, wtilde_array, my_igp, n1, aqsmtemp, aqsntemp, I_eps_array, ssxt, scht, ncouls, igp, scha);

            sch_array[iw] +=GPPComplex_mult(scht, 0.5);
        }

        for(int iw=nstart; iw<nend; ++iw)
            achtemp[iw] += GPPComplex_mult(sch_array[iw] , vcoul[igp]);//Store final output here

        acht_n1_loc[n1] += GPPComplex_mult(sch_array[2] , vcoul[igp]);
    }
}

int main(int argc, char** argv)
{

//The input to the executable needs 4 arguments.
    if (argc != 5)
    {
        std::cout << "The correct form of input is : " << endl;
        std::cout << " ./a.out <number_bands> <number_valence_bands> <number_plane_waves> <matrix_divider> " << endl;
        exit (0);
    }

//Input parameters stored in these variables.
    const int number_bands = atoi(argv[1]);
    const int nvband = atoi(argv[2]);
    const int ncouls = atoi(argv[3]);
    const int nodes_per_group = atoi(argv[4]);

//Constants that will be used later
    const int npes = 1; 
    const int ngpown = ncouls / (nodes_per_group * npes); 
    const double e_lk = 10;
    const double dw = 1;
    const double to1 = 1e-6;
    const double gamma = 0.5;
    const double sexcut = 4.0;
    const double limitone = 1.0/(to1*4.0);
    const double limittwo = pow(0.5,2);
    const double e_n1kq= 6.0; 
    const double occ=1.0;

    //Printing out the params passed.
    std::cout << "number_bands = " << number_bands \
        << "\t nvband = " << nvband \
        << "\t ncouls = " << ncouls \
        << "\t nodes_per_group  = " << nodes_per_group \
        << "\t ngpown = " << ngpown \
        << "\t nend = " << nend \
        << "\t nstart = " << nstart \
        << "\t gamma = " << gamma \
        << "\t sexcut = " << sexcut \
        << "\t limitone = " << limitone \
        << "\t limittwo = " << limittwo << endl;

    // Memory allocation of input data structures.
    // Two dimensional arrays from theory have been initialized as a single dimension in m*n format for performance.
    GPPComplex *acht_n1_loc = new GPPComplex [number_bands];
    GPPComplex *aqsmtemp = new GPPComplex [number_bands*ncouls];
    GPPComplex *aqsntemp = new GPPComplex [number_bands*ncouls];
    GPPComplex *I_eps_array = new GPPComplex [ngpown*ncouls];
    GPPComplex *wtilde_array = new GPPComplex [ngpown*ncouls];
    int *inv_igp_index = new int[ngpown];
    double *vcoul = new double[ncouls];
    double wx_array[nend-nstart];

    //arrays that will be later used to store the output results
    GPPComplex achtemp[nend-nstart]; 
    GPPComplex asxtemp[nend-nstart];

    GPPComplex achstemp = GPPComplex(0.0, 0.0);

    //Data structures that store intermediete results
    GPPComplex ssx_array[nend-nstart], \
        sch_array[nend-nstart], \
        scht, ssxt;
    GPPComplex *ssxa = new GPPComplex [ncouls];
    GPPComplex *scha = new GPPComplex [ncouls];

    //Printing the size of each of the input data structures.
    cout << "Size of wtilde_array = " << (ncouls*ngpown*2.0*8) / pow(1024,2) << " Mbytes" << endl;
    cout << "Size of aqsntemp = " << (ncouls*number_bands*2.0*8) / pow(1024,2) << " Mbytes" << endl;
    cout << "Size of I_eps_array array = " << (ncouls*ngpown*2.0*8) / pow(1024,2) << " Mbytes" << endl;

    //Some expressions declared to be used later in the initialization.
    GPPComplex expr0( 0.0 , 0.0);
    GPPComplex expr( 0.5 , 0.5);

//Initializing the data structures
   for(int i=0; i<number_bands; i++)
       for(int j=0; j<ncouls; j++)
       {
           aqsntemp[i*ncouls+j] = GPPComplex_mult(expr, (double)(i+j));
           aqsmtemp[i*ncouls+j] = GPPComplex_mult(expr, (double)(i+j));
       }


   for(int i=0; i<ngpown; i++)
   {
       for(int j=0; j<ncouls; j++)
       {
           I_eps_array[i*ncouls+j] = GPPComplex_mult(expr, (double)(i+j));
           wtilde_array[i*ncouls+j] = GPPComplex_mult(expr, (double)(i+j));
       }

        inv_igp_index[i] = (i+1) * ncouls / ngpown;
   }

   for(int i=0; i<ncouls; i++)
   {
       ssxa[i] = expr0;
       scha[i] = expr0;
       vcoul[i] = 1.0*i;
   }


    for(int iw=nstart; iw<nend; ++iw)
    {
        achtemp[iw] = expr0;
        asxtemp[iw] = expr0;
        wx_array[iw] = e_lk - e_n1kq + dw*((iw+1)-2);
        if(wx_array[iw] < to1) wx_array[iw] = to1;
    }

    //Start the timer before the work begins.
    timeval startTimer, endTimer;
    gettimeofday(&startTimer, NULL);


    for(int n1 = 0; n1<nvband; ++n1) 
    {
        reduce_achstemp(n1, inv_igp_index, ncouls, aqsmtemp, aqsntemp, I_eps_array, achstemp, ngpown, vcoul);

        for(int my_igp=0; my_igp<ngpown; ++my_igp)
        {
            int igp = inv_igp_index[my_igp];
            if(igp >= ncouls)
                igp = ncouls-1;

            //Reinitialize the intermediete variables to initial values.
            for(int i=nstart; i<nend; i++)
            {
                ssx_array[i] = expr0;
                sch_array[i] = expr0;
            }

            for(int iw=nstart; iw<nend; iw++)
            {
                scht = ssxt = expr0;
                double wxt = wx_array[iw];
                flagOCC_solver(wxt, wtilde_array, my_igp, n1, aqsmtemp, aqsntemp, I_eps_array, ssxt, scht, ncouls, igp, ssxa, scha);
                ssx_array[iw] += ssxt;
                sch_array[iw] += GPPComplex_mult(scht, 0.5);
                
                asxtemp[iw] += GPPComplex_mult(ssx_array[iw] , occ * vcoul[igp]);//Store output of the first nvband iterations.
            }

            for(int iw=nstart; iw<nend; ++iw)
                achtemp[iw] += GPPComplex_mult(sch_array[iw] , vcoul[igp]);//Store final output here

            acht_n1_loc[n1] += GPPComplex_mult(sch_array[2] , vcoul[igp]);
        }
    }

    #pragma acc data copyin(vcoul[0:ncouls])
    #pragma acc parallel loop
    for(int n1 = nvband; n1<number_bands; ++n1) 
    {
        reduce_achstemp(n1, inv_igp_index, ncouls, aqsmtemp, aqsntemp, I_eps_array, achstemp, ngpown, vcoul);

        for(int my_igp=0; my_igp<ngpown; ++my_igp)
        {
            int igp = inv_igp_index[my_igp];
            if(igp >= ncouls)
                igp = ncouls-1;

            //Reinitialize the intermediete variables to initial values.
            for(int i=nstart; i<nend; i++)
            {
                ssx_array[i] = expr0;
                sch_array[i] = expr0;
            }

            for(int iw=nstart; iw<nend; ++iw)
            {
                    scht = ssxt = expr0;
                    double wxt = wx_array[iw];
                    noflagOCC_solver(wxt, wtilde_array, my_igp, n1, aqsmtemp, aqsntemp, I_eps_array, ssxt, scht, ncouls, igp, scha);

                    sch_array[iw] += GPPComplex_mult(scht, 0.5);
            }
            for(int iw=nstart; iw<nend; ++iw)
                achtemp[iw] += GPPComplex_mult(sch_array[iw] , vcoul[igp]);//Store final output here

            acht_n1_loc[n1] += GPPComplex_mult(sch_array[2] , vcoul[igp]);
        }
    }


    //Time Taken
    gettimeofday(&endTimer, NULL);
    double elapsedTimer = (endTimer.tv_sec - startTimer.tv_sec) +1e-6*(endTimer.tv_usec - startTimer.tv_usec);


    for(int iw=nstart; iw<nend; ++iw)
        achtemp[iw].print();

    cout << "********** Time Taken **********= " << elapsedTimer << " secs" << endl;

    //Free the allocated memory
    free(acht_n1_loc);
    free(wtilde_array);
    free(aqsmtemp);
    free(aqsntemp);
    free(I_eps_array);
    free(inv_igp_index);
    free(vcoul);

    return 0;
}
