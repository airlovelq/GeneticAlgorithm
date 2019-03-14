#ifndef _GENETICALGORITHM_H
#define _GENETICALGORITHM_H


#include <stdio.h>

//染色体结构体
struct Chromosome
{
	int *bit;           //码串
	int nlen;             //码串长度
	double fitness;       //适应值
	double rate_fit;      //相对的fit值，即所占百分比，用于轮盘赌算法
	double cumu_fit;
	Chromosome()
	{
		bit = NULL;
		fitness = 0.0;
		rate_fit = 0.0;
		cumu_fit = 0.0;
		nlen = 0;
	}
	~Chromosome()
	{
		if (bit)
			delete bit;
	}
	Chromosome(const Chromosome &a) throw()
	{
		fitness = a.fitness;
		rate_fit = a.rate_fit;
		cumu_fit = a.cumu_fit;
		nlen = a.nlen;
		//if (bit)
			//delete bit;
		bit = new int[nlen];
		for ( int i=0 ; i<nlen ; i++ )
		{
			bit[i] = a.bit[i];
		}
	}
	Chromosome& operator= ( const Chromosome& a) throw()
	{
		fitness = a.fitness;
		rate_fit = a.rate_fit;
		cumu_fit = a.cumu_fit;
		nlen = a.nlen;
		if (bit)
			delete bit;
		bit = new int[nlen];
		for ( int i=0 ; i<nlen ; i++ )
		{
			bit[i] = a.bit[i];
		}
		return *this;
	}
	void clone(Chromosome &n) throw()
	{
		n.fitness = fitness;
		n.cumu_fit = cumu_fit;
		n.rate_fit = rate_fit;
		n.nlen = nlen;
		if ( n.bit )
			delete n.bit;
		n.bit = new int[nlen];
		for ( int i=0 ; i<nlen ; i++ )
		{
			n.bit[i] = bit[i];
		}
	}
	void Init(int len) throw()
	{
		nlen = len;
		bit = new int[len];
	}

};

class  CGeneticAlgorithmBase
{
public:
	CGeneticAlgorithmBase(int *poprange,int nrangelen,double dcross_rate=0.9,double dmutate_rate=0.05,int npopulation=500,int niterations=100,int nDNA_size=5,int nTerminateSec = 20 );
	virtual ~CGeneticAlgorithmBase(void);

protected:

	double m_dcross_rate;    //交配概率
	double m_dmutate_rate;   //基因突变概率
	int m_npopulation;    //种群大小
	int m_niterations;    //迭代次数
	int m_nDNA_size;      //DNA长度

	Chromosome *m_spop;           //样本

	int *m_poprange;            //从该集合中抽取样本
	int m_nrangelen;               //集合长度

	int m_nTerminateSec;      //终止时间
	Chromosome m_best;            //最优个体

protected:
	

	virtual bool InitPopulation();         //生成样本
	virtual bool Objective() = 0;            //计算适应度值
	virtual bool Cross();                 //交叉操作
	virtual bool Mutation();              //变异操作
	virtual bool seletc_prw();            //轮盘赌选择操作
	virtual bool fresh_property();        //函数：更新种群内个体的属性值

public:
	virtual bool Generate();          //执行过程

	Chromosome getBest() const;       //获取最优个体
private:
	CGeneticAlgorithmBase(){};
	CGeneticAlgorithmBase(const CGeneticAlgorithmBase &a){};
    CGeneticAlgorithmBase& operator = (const CGeneticAlgorithmBase &a){return *this;};

};

#endif