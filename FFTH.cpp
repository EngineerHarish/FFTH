#include <iostream>
#include <cmath>
#include <complex>
#include <iomanip>

using namespace std;
//Complex datatype
typedef complex<float> cox;

//Declare functions
void getSigcount(int &N);
void getSig(cox *sig, int N);
void getOp(int &mode);
void getFlow();
void runAlgorithm(cox &sig, int &mode, int N);
void algoDIF(cox &sig, int N, int mode);
void algoDIT(cox &sig, int N, int mode);
void inverseOrnot(cox *sig, int N);
void bitReverse(cox *sig, int N);
void butterfly(cox *sig, int N, int mode);
void printSig(cox *sig, int N, int j);
void twiddle (cox *twid, int M);
void reset (int &index);


//Global variables
int choice; //Stores the choice of the user for FFT or inverse FFT

//Main function
int main()
{   
    int N; //Number of inputs
    int mode; //Stores the algorithm to be performed

    //Functions
    getSigcount(N);//Get number of inputs
    
    //Declare the complex signal array
    cox* sig= new cox[N];

    getFlow();//Get if the user wants to perform the inverse FFT or just the FFT

    getSig(sig, N);//Receives the input data

    getOp(mode);//Get the algorithm to be performed

    runAlgorithm(*sig, mode, N);//Run the algorithm chosen by the user

    delete[]sig;
    return 0;
}

//Number of inputs function
void getSigcount(int &N)
{
    cout << "Enter the number of inputs: ";
    cin>> N;
}

//Input array function
void getSig(cox *sig, int N)
{   
    float realPart, imagPart;//tmp Variables to receive input terms

    for (int i= 0; i< N; i++)
    {
        //Get the real value
        cout<< "Enter the real term of the element " << i+1 << ": ";
        cin >> realPart;

        //Get the imaginary value
        cout<< "Enter the imaginary term of the element " << i+1 << ": ";
        cin >> imagPart;

        //Assign that value to the input array
        sig[i]= cox(realPart, imagPart);
        cout << endl;
    }
}

void getFlow()
{
    //Get if the user wants to perform the inverse FFT or just the FFT
    cout<< "Do you want to perform FFT or inverse FFT? (1 for FFT, 0 for inverse FFT): ";
    cin >> choice;
}

//Get the algorithm to be performed
void getOp(int &mode)
{   
    //List out the algorithm options
    cout << "********************Enter the algorithm to be used********************" <<endl;
    cout << "1) DIF\n2)DIT" <<endl;

    cout << "Select the algorithm you want to perform: ";

    //Get the algorithm to be used
    cin >> mode;
}

//Algorithm selection function
void runAlgorithm(cox &sig, int &mode, int N)
{
    switch(mode)
    {
        case 1: 
            algoDIF(sig, N, mode);
        break;
        
        case 2:
            algoDIT(sig, N, mode);
        break;
    }  
}

//DIF algorithm fucntion
void algoDIF(cox &sig, int N, int mode)
{    
    //Butterfly operation
    butterfly(&sig, N, mode);

    //Bit reversal
    bitReverse(&sig, N);
}

//DIF algorithm fucntion
void algoDIT(cox &sig, int N, int mode)
{
    //Bit reversal
    bitReverse(&sig, N);

    //Butterfly operation
    butterfly(&sig, N, mode);
}

//Bit Reversal function
void bitReverse(cox *sig, int N)
{   
    int BW= log2(N);

    cox* revSig= new cox[N];

    //Bit reversal operation
    for (int i = 0; i < N; ++i) 
    {
        int revIndex = 0;
        for (int j = 0; j < BW; ++j)
        {
            revIndex |= ((i >> j) & 1) << (BW - 1 - j);
        }
        if (revIndex < N) 
        {
            revSig[revIndex] = sig[i];
        }
    }

    for (int i= 0; i < N; i++)
    {
        sig[i]= revSig[i];
    }
}

//Butterfly function
void butterfly(cox *sig, int N, int mode)
{   

    int flipper;
    int reach;
    int check;

    cox* tmp= new cox[N];
    //Variables for Twiddle Factor
    cox *twid= new cox[N];

    
    switch (mode)
    {   
        //Buttefly operation for DIF algorithm
        case 1:

            flipper= 1;
            reach= N/2;

            

            for (int j= log2(N); j >= 1; j--)
            {   
                //Twiddle Factor
                twiddle (twid, N);
                int index= 0;

                //Buterfly operation implementation
                for (int i= 0; i<N; i++)
                {
                    if(flipper== 1 )
                    {
                        tmp[i]= sig[i]+sig[i+reach];

                        //Check if the index has reached the end of the current reach
                        check = (i+1) % reach;

                        //If it has, flip the flipper
                        if(check == 0)
                        {
                            flipper*= -1;
                        }
                        
                    }
                    else
                    {
                        tmp[i]= (-sig[i]+sig[i-reach])*twid[index];
                        
                        index++; //Increment the Twiddle index

                        //Check if the index has reached the end of the current reach
                        check = (i+1) % reach;

                        //If it has, flip the flipper and reset the Twiddle index
                        if(check == 0)
                        {
                            flipper*= -1;

                            reset (index);
                        }
                    }

                    if (i== N-1)
                    {   
                        for (int k= 0; k<N; k++)
                            {
                                sig[k]= tmp[k];
                            }

                        reach /= 2;
                    }
                
                    
                }

                //Special case for the last stage. Array is reversed for the DIF algorithm.
                if (j == 1)
                {
                    //Reverse the array
                    bitReverse(sig, N);

                    if (choice==0)
                    {
                        //If the user wants to perform the inverse FFT, divide the signal by N
                        inverseOrnot(sig, N);
                    }

                    //Then print it
                    printSig(sig, N, j);
                }
                    
                else
                    printSig(sig, N, j);
            }
            
        break;

        //Butterfly diagram for DIT algorithm
        case 2:

            flipper= 1;
            reach= 1;

            //Full butterfly logic with a time complexity of O(N log(N))
            for (int j= 1; j <= log2(N); j++)
            {   
                //Twiddle Factor
                twiddle (twid, N);
                int index= 0;

                //Buterfly operation implementation
                for (int i= 0; i<N; i++)
                {
                    if(flipper== 1)
                    {
                        tmp[i]= sig[i]+(sig[i+reach]*twid[index]);

                        index++; //Increment the Twiddle Index

                        check = (i+1) % reach;
                        if(check == 0)
                        {
                            flipper*= -1;

                            reset (index);
                        }
                        
                    }
                    else
                    {
                        tmp[i]= (-sig[i]*twid[index])+sig[i-reach];

                        index++; //Increment the Twiddle index
                        
                        check = (i+1) % reach;
                        if(check == 0)
                        {
                            flipper*= -1;

                            reset (index);
                        }
                    }

                    if (i== N-1)
                    {   
                        for (int k= 0; k<N; k++)
                            {
                                sig[k]= tmp[k];
                            }

                        reach *= 2;
                    }
                }

                if (choice==0 && j == log2(N))
                {
                    //If the user wants to perform the inverse FFT, divide the signal by N
                    inverseOrnot(sig, N);
                }
                
                //Print output for the current stage
                printSig(sig, N, j);
            }

        break;
    
    delete [] twid;
    delete []tmp; 
    }
}

//Twiddle Factor Function
void twiddle (cox *twid, int M)
{
    float Pi= 3.1415;

    int n= log2(M);
    
    for (int k= 0; k<n; k++)
    {

        int newN= pow(2,k);

        int sub= pow(2, k+1);

        for (int j=0; j<newN; j++)
        {
            float theta = -2.0 * Pi * j/sub;

            twid[j]= cox (cos(theta), sin(theta));
        }
    }
}

//Twiddle Index Reset function
void reset (int &index)
{
    index= 0;
}

void inverseOrnot(cox *sig, int N)
{
        for (int i= 0; i<N; i++)
        {
            sig[i]/= N; //Divide the signal by N to get the inverse FFT
        }
}

//Print the signal for the current Stage
void printSig(cox *sig, int N, int j)
{
    cout << "\nOutput of stage " << j <<endl;
                
    //Print the butterflied output (Increases time complexity)
    for (int i= 0; i<N; i++)
    {   
        //You can use this same loop at the end of all of the stages
        cout << setprecision(4) << sig[i] << endl ; //Set precision is used for better readability
    }  
}