// standard include
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <windows.h>
using namespace std;

// define priority classes
#define NORMAL_PRIORITY_CLASS       0x00000020
#define IDLE_PRIORITY_CLASS         0x00000040
#define HIGH_PRIORITY_CLASS         0x00000080
#define REALTIME_PRIORITY_CLASS     0x00000100

// define parameters
#define L           100      /* lattice size                   */
// #define SIZE        (L*L)    /* number of sites                */
#define SIZE        100    /* number of sites                */
#define MC_STEPS    61000   /* run-time in MCS     */
#define K           0.1      /* temperature */
#define NAMEOUT     "K4b075r5Q2"
#define RANDOMIZE   3145215

int defector1,cooperator1,defector2,cooperator2;
double score;
double b;
double delta;
const double Upper=1.2;
const double Lower=0.0;


typedef int       tomb1[SIZE];
typedef long int  tomb3[SIZE][SIZE];
typedef double    tomb4[SIZE];
typedef long double  tomb5[SIZE][SIZE];

tomb1 player_s;           /* matrix ,containing players strategies */
tomb1 player_a;           /* matrix ,containing players strategies */
tomb3 player_n;           /* matrix, containing players neighbours */
tomb5 player_we;           /* matrix, containing players contribute */
tomb3 player_adjacency;     //邻接矩阵
tomb1 player_degree;        //每个节点的度
tomb4 player_payoff;          //收益矩阵

void prodgraph(void);      // 初始化每个个体的邻居,即把每个个体邻居的编号保存起来
void initial(void);        // 初始化空间中每个个体的策略
void update(void);          // 策略更新
void tongji(void);          // 统计群体中策略的数目
double payoff(int i);       // 个体i收益计算
int degree(int i);          // 计算节点 i 的度
void payoffCalculate(void); // 计算收益，放入矩阵player_payoff中
void degreeCalculate(void); // 计算节点度，放入player_degree中
void createLink(void);      //建立新的连接

ofstream outfile1;
ofstream outfile2;
ofstream outfile3;

/*************************** RNG procedures ****************************************/
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[NN]; /* the array for the state vector  */
static int mti=NN+1; /* mti==NN+1 means mt[NN] is not initialized */
void sgenrand(unsigned long seed)
{int i;
 for (i=0;i<NN;i++) {mt[i] = seed & 0xffff0000; seed = 69069 * seed + 1;
                     mt[i] |= (seed & 0xffff0000) >> 16; seed = 69069 * seed + 1;
  }
  mti = NN;
}
void lsgenrand(unsigned long seed_array[])
{ int i; for (i=0;i<NN;i++) mt[i] = seed_array[i]; mti=NN; }
double genrand()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    if (mti >= NN)
    {
        int kk;
        if (mti == NN+1) sgenrand(4357);
        for (kk=0;kk<NN-MM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MM] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<NN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MM-NN)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[NN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[NN-1] = mt[MM-1] ^ (y >> 1) ^ mag01[y & 0x1];
        mti = 0;
    }
    y = mt[mti++]; y ^= TEMPERING_SHIFT_U(y); y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C; y ^= TEMPERING_SHIFT_L(y);
    return y;
}

double randf(){ return ( (double)genrand() * 2.3283064370807974e-10 ); }	//随机浮点数
long randi(unsigned long LIM){ return((unsigned long)genrand() % LIM); }	//随机整数

/********************** END of RNG ************************************/

/*********************初始化空间个体中的策略**************************/
void initial(void)
{
	 int i,j;
	 int ii;
	 int temp;
	 int iii,jjj;
     cooperator1=0;
	 defector1=0;
    for (i=0; i<SIZE; i++)
    {
        player_s[i]=(int)randi(2);
        switch(player_s[i])
        {
            case 0: cooperator1++;  break;
            case 1: defector1++;    break;
        }

    }
}

/*************初始化邻接矩阵**************/
void initialAdjacency(int proportion)
{
    int i,j;
    int player1, player2;
    for(i=0; i<SIZE; i++)                   //需要先清零
    {
        for(j=0; j<SIZE; j++)
        {
            player_adjacency[i][j]=0;
        }
    }

    for (i = 0; i<SIZE ;i++)                //按照比例分配是否连接
    {
        player1 = i;
        for (j=0; j<SIZE; j++)
        {
            if(i != j)
            {
                if (randi(100)<proportion)
                    player_adjacency[player1][j]=1;
                else
                    player_adjacency[player1][j]=0;
            }

        }
    }
}

/************************* 个体收益计算 **************************/
double payoff(int i)
{
	int star1,star2;
	int player1,player2;
    double score=0.0;
    double score1=0.0;
	int j;
	player1=i;
	star1=player_s[player1];
		for(j=0; j<SIZE; j++)
		{
            if (j != player1)
            {
                player2 = j;
                star2 = player_s[player2];
                switch(2*star1+star2)
                {
                    case 0:score=score+player_adjacency[player1][j]*1.0;break;
                    case 1:score=score+player_adjacency[player1][j]*0.0;break;
                    case 2:score=score+player_adjacency[player1][j]*b;break;
                    case 3:score=score+player_adjacency[player1][j]*0.0;break;
                }
            }
            score1 += score;
		}
	return score1;
}

/************** 计算每个节点的收益 ******************/
void payoffCalculate(void)
{
    int j;
    for( j=0; j<SIZE; j++)
    {
        player_payoff[j] = payoff(j); 
    }
}

/*************** 计算该节点i的度 *********************/
int degree(int i)
{
    int player1;
    int k , j;
    k = 0;
    player1 = i;
    for (j=0; j<SIZE; j++)
    {
        if (player_adjacency[player1][j] != 0)
        {
            k++;
        }
    }
    return k;
}

void degreeCalculate(void)
{
    int i;
    for(i=0; i<SIZE; i++)
    {
        player_degree[i] = degree(i);
    }
}

/****************建立新连接*********************/
void createLink(void)
{
    int i;
    int player1, player2;
    double femi=0.0;
    double sumGain=0.0;
    double averageGain=0.0;
    do
    {
        player1 = (int)randi(SIZE);
        do
        {
            player2 = (int)randi(SIZE);
        }while(player1 == player2);
    }while( player_adjacency[player1][player2] != 0 );

    for(i=0; i<SIZE; i++)
    {
        sumGain += player_payoff[i];
    }
    averageGain = sumGain/SIZE;

    femi = abs(player_payoff[player1]-player_payoff[player2]) / (2*averageGain);
    //if(averageGain == 0.0) femi = 1;                                            //避免最开始出现0/0的情况
    if (femi > randf())
    player_adjacency[player1][player2] = 1;
    player_adjacency[player2][player1] = 1;

}

/****************************** 策略更新 **************************/
void update(void)
{
	int i,j,k,choice;
    int player1,player2,player3;
    int degree, degree1, degree2;
    int numberOfNeighbour=0;                            //统计邻居有多少个
    int numberOfCount=0;                                //统计数到了第几个邻居
    double numberOfRandom=0.0;                          //生成的随机数
	double score1,score2,dscore;
	double femi;
	for(i=0;i<SIZE;i++)
	{
        numberOfCount = 0;
        numberOfNeighbour = 0;
        numberOfRandom = 0.0;
		player1 = (int) randi(SIZE);					//为什么是随机选的player1？？？？？？？？？？？？？
		// Weight_update(player1);							//异步更新，所以随机选择，然后平均每轮没人有一次机会就可以
        
        for(j=0; j<SIZE; j++)                           //统计有多少个邻居
        {
            if(player_adjacency[player1][j] != 0)
            {
                numberOfNeighbour++;
            }
        }
        if (numberOfNeighbour != 0)
        {
            numberOfRandom = randi(numberOfNeighbour) + 1;  //随机挑选邻居
            for(k=0; k<SIZE; k++)
            {
                if (player_adjacency[player1][k] != 0)
                {
                    numberOfCount++;
                    if (numberOfCount >= numberOfRandom)
                    {
                        choice = k;break;                   //挑选好之后推出循环
                    }
                }
            }
            player2=choice;
    
            degree1 = player_degree[player1];              //找到度大的节点
            degree2 = player_degree[player2]; 
            if (degree1 > degree2)
                degree = degree1;
            else
                degree = degree2;
            score1 = player_payoff[player1];
            score2 = player_payoff[player2];
            if(player_s[player1]!=player_s[player2])
            {
                dscore=score2/degree-score1/degree;				//在这一步中已经确定了，为负时候肯定不变，为正再看概率
                femi=dscore/b;
                if (femi>randf()) player_s[player1]=player_s[player2];
            }
        }
          
	}

}

/************************* 统计群体中各种策略的数目 **************************/
void tongji(void)
{
	int i;
	cooperator1=0;
	defector1=0;
	for(i=0;i<SIZE;i++)
	{
		switch(player_s[i])
		{
            case 0: cooperator1++; break;
            case 1: defector1++; break;
		}
	}
}


/******************************  主函数   ***************************/
int main()
{
    int steps, propor;
    int countNumber = 0;
    int i,j;
	double x,XX,aa,x1,XX1,aa1;
	outfile1.open("frequency.txt");
    outfile2.open("average.txt");
    outfile3.open("degree.csv");
	if(!outfile1)
	{
		cout<<"can not open";
		abort();
	}
	if(!outfile2)
	{
		cout<<"can not open";
		abort();
    }
    if(!outfile3)
    {
        cout<<"can not open degree.csv";
        abort();
    }

	// initialize the random number generation
	sgenrand(RANDOMIZE);

    cout<<"Start the experiment!"<<endl;
	// begins the mc steps
	int ii;    
    for(b=1.0; b<1.8 ; b=b+0.02)
    {
        initial();                                      //初始化空间个体中的策略
        for ( propor=10; propor<60; propor=propor+10)    //分别从0-50%的初始状态进行迭代
        {
            cout<<"The temptation to defect is :"<<b<<" and the proportion of beginning is:"<<propor<<endl;
            aa=0;
            aa1=0;
            initialAdjacency(propor);
            for(steps=0; steps<MC_STEPS; steps++)
            {
                
                payoffCalculate();                      //计算节点收益
                degreeCalculate();                      //计算节点的度
                createLink();                           //添加新连接
                update();                               //更新策略
                tongji();
                x = (double)cooperator1/SIZE;
                x1 = (double)defector1/SIZE;            
                outfile1<<steps<<'\t'<< x <<'\t'<< x1 <<endl;//C,D
                if((x==1)||(x1==1))
                {
                    XX = x;
                    XX1 = x1;
                    //break;                             
                }

                for(i=0; i<SIZE; i++)                   //矩阵满时候停止本轮循环
                {
                    for(j=0; j<SIZE; j++)
                    {
                        if (player_adjacency[i][j] !=0)
                            countNumber++;
                    }
                }
                if (countNumber == SIZE*SIZE)
                    break;

                if(steps > MC_STEPS-1001)
                {
                    aa += x;
                    aa1 += x1;
                    XX = (double)(aa)/1000;
                    XX1 = (double)(aa1)/1000;
                }
            }
            outfile2<<propor<<'\t'<<b<<'\t'<<XX<<'\t'<<XX1<<endl;//C，D
            cout<<propor<<'\t'<<b<<'\t'<<XX<<'\t'<<XX1<<endl;//最后输出待定?????????????????????????

            outfile3<<propor<<','<<b<<endl;
            for(int num=0; num<SIZE; num++)
            {
                outfile3<<player_degree[num]<<',';
                // cout << player_degree[num]<<'\t';
            }
            outfile3<<endl;
            // cout<<endl;
        }
    }
	outfile1.close();
    outfile2.close();
    outfile3.close();
	printf("This is the end\n");
	return 0;
}