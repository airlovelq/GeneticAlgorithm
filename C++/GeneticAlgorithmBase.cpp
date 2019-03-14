#include "GeneticAlgorithmBase.h"
#include <random>
#include <time.h>

CGeneticAlgorithmBase::CGeneticAlgorithmBase(int *poprange,int nrangelen,double dcross_rate, double dmutate_rate, int npopulation, int niterations, int nDNA_size,int nTerminateSec )
    :m_poprange(poprange),m_nrangelen(nrangelen),m_dcross_rate(dcross_rate),m_dmutate_rate(dmutate_rate),m_npopulation(npopulation),m_niterations(niterations),m_nDNA_size(nDNA_size),m_nTerminateSec(nTerminateSec)
{
	m_poprange = new (std::nothrow) int[nrangelen];
	memcpy(m_poprange,poprange,nrangelen*sizeof(int));
}


CGeneticAlgorithmBase::~CGeneticAlgorithmBase()
{
	if ( m_poprange )
		delete []m_poprange;
}

bool CGeneticAlgorithmBase::InitPopulation() 
{

	m_spop = new (std::nothrow) Chromosome[m_npopulation];
	if ( m_spop == NULL )
		return false;
	srand((unsigned)time(NULL));
	for ( int i = 0; i < m_npopulation; i++ )
	{
		m_spop[i].Init(m_nDNA_size);
		for (int j = 0; j < m_nDNA_size; j++)
		{
			int xh = rand()%(m_nrangelen);
			m_spop[i].bit[j] = m_poprange[xh];
			m_spop[i].fitness = 1.0;
		}
	}
	return true;
}


bool CGeneticAlgorithmBase::Cross()
{
	int i = 0, j = 0;
	double rate_rand = 0.0;
	short int bit_temp = 0;
	int num1_rand = 0, num2_rand = 0, position_rand = 0;
	//产生随机数种子
	srand((unsigned)time(NULL));
	//应当交叉变异多少次呢？先设定为种群数量
	for (j = 0; j<m_npopulation; j++)
	{
		rate_rand = (double)rand() / (RAND_MAX);
		//如果大于交叉概率就进行交叉操作
		if (rate_rand <= m_dcross_rate)
		{
			//随机产生两个染色体
			num1_rand = (int)rand() % (m_npopulation);
			num2_rand = (int)rand() % (m_npopulation);
			//随机产生两个染色体的交叉位置
			position_rand = (int)rand() % (m_nDNA_size);
			//采用单点交叉，交叉点之后的位数交换
			for (i = position_rand; i<m_nDNA_size; i++)
			{
				bit_temp = m_spop[num1_rand].bit[i];
				m_spop[num1_rand].bit[i] = m_spop[num2_rand].bit[i];
				m_spop[num2_rand].bit[i] = bit_temp;
			}

		}
	}
	return true;
}

// 函数：变异操作
bool CGeneticAlgorithmBase::Mutation()
{
	int position_rand = 0;
	int i = 0;
	double rate_rand = 0.0;
	//产生随机数种子
	srand((unsigned)time(NULL));
	//变异次数设定为种群数量
	for (i = 0; i<m_npopulation; i++)
	{
		rate_rand = (float)rand() / (RAND_MAX);
		//如果大于交叉概率就进行变异操作
		if (rate_rand <= m_dmutate_rate)
		{
			//随机产生突变位置
			position_rand = (int)rand() % (m_nDNA_size);
			//double dmutationvalue = m_poprange[int(double(rand()) / RAND_MAX*(m_nrangelen - 1))];
			int xh = rand()%(m_nrangelen);
			//突变
			//m_spop[i].bit[position_rand] = !m_spop[i].bit[position_rand];
			m_spop[i].bit[position_rand] = m_poprange[xh];
		}

	}
	return true;
}

bool CGeneticAlgorithmBase::seletc_prw()
{
	Chromosome *spopcur = new (std::nothrow) Chromosome[m_npopulation];
	if ( !spopcur )
		return false;
	for (int k = 0; k < m_npopulation; k++)
	{
		m_spop[k].clone(spopcur[k]);
	}
	
	int i = 0, j = 0;
	double rate_rand = 0.0;
	m_best = spopcur[0];
	//产生随机数种子
	srand((unsigned)time(NULL));
	for (i = 0; i < m_npopulation; i++)
	{
		rate_rand = (float)rand() / (RAND_MAX);
		if (rate_rand < spopcur[0].cumu_fit)
		{
			m_spop[i] = spopcur[0];
			m_spop[i].fitness = 0.0;
			m_spop[i].cumu_fit = 0.0;
			m_spop[i].rate_fit = 0.0;
		}
		else
		{
			for (j = 0; j < m_npopulation-1; j++)
			{
				if (spopcur[j].cumu_fit <= rate_rand && rate_rand < spopcur[j + 1].cumu_fit)
				{
					m_spop[i] = spopcur[j+1];
					m_spop[i].fitness = 0.0;
					m_spop[i].cumu_fit = 0.0;
					m_spop[i].rate_fit = 0.0;
					break;
				}
			}
		}

		//如果当前个体比目前的最有个体还要优秀，则将当前个体设为最优个体
		if (spopcur[i].fitness > m_best.fitness)
			m_best = spopcur[i];
	}
	delete []spopcur;
	return true;

}

bool CGeneticAlgorithmBase::fresh_property()
{
	int j = 0;
	double sum = 0.0;
	bool bRes = Objective();
	if ( !bRes )
		return bRes;
	for (j = 0; j < m_npopulation; j++)
	{

		//染色体解码，将二进制换算为十进制，得到一个整数值
		//计算二进制串对应的10进制数值
		//decode(population_current[j]);
		//计算染色体的适应度
		//population_current[j].fitness = Objective(m_spop[j].bit);
		sum = sum + m_spop[j].fitness;

	}


	//计算每条染色体的适应值百分比及累计适应度值的百分比，在轮盘赌选择法时用它选择染色体  
	m_spop[0].rate_fit = m_spop[0].fitness / sum;
	m_spop[0].cumu_fit = m_spop[0].rate_fit;
	for (j = 1; j < m_npopulation; j++)
	{
		m_spop[j].rate_fit = m_spop[j].fitness / sum;
		m_spop[j].cumu_fit = m_spop[j].rate_fit + m_spop[j - 1].cumu_fit;
	}
	return bRes;

}

bool CGeneticAlgorithmBase::Generate()
{
	bool bRes = InitPopulation();
	if ( !bRes )
		return bRes;
	bRes = fresh_property();
	if ( !bRes )
		return bRes;
	clock_t tstart = clock();
	for (int i = 0; i < m_niterations; i++)
	{
		clock_t tend = clock();
		if ( double(tend-tstart)/1000.0 >= m_nTerminateSec )
			break;
		//挑选优秀个体组成新的种群
		bRes = seletc_prw();
		if ( !bRes )
			break;
		//对选择后的种群进行交叉操作
		bRes = Cross();
		if ( !bRes )
			break;
		//对交叉后的种群进行变异操作
		bRes = Mutation();
		if ( !bRes )
			break;
		//更新种群内个体的属性值
		bRes = fresh_property();
		if ( !bRes )
			break;
	}
	return bRes;

}

Chromosome CGeneticAlgorithmBase::getBest() const
{
	return m_best;
}

