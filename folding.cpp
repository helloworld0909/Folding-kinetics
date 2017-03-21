#include "stdafx.h"
#include "Struct_Chain.h"
#define NUMTEST 100
#define τ 100000
#define tINACF 1000000
#define ACFSTEP 1000

int InitNode(Node chain[], int n, int typeArray[]) {//产生笔直向上的链
	if (chain == NULL) {
		cout << "Allocate memory fail" << endl;
		exit(1);
	}
	for (int i = 1; i < n - 1; i++) {
		chain[i].type = typeArray[i];
		chain[i].nextPos = 1;
		chain[i].prevPos = -1;
	}
	chain[0].type = typeArray[0];
	chain[0].prevPos = 0;
	chain[0].nextPos = 1;
	chain[n - 1].type = typeArray[n - 1];
	chain[n - 1].nextPos = 0;
	chain[n - 1].prevPos = -1;
	/*if (checkchain(chain)) {
		cout << "Init Error" << endl;
		print(chain);
		exit(1);
	}*/
	return 1;
}
int InitNode(Node chain[], int n, int PosArray[], int typeArray[]) {//输入的PosArray有n-1个，对应nextPos
	if (chain == NULL) {
		cout << "Allocate memory fail" << endl;
		exit(1);
	}
	for (int i = 0; i < n - 1; i++) {
		chain[i].type = typeArray[i];
		chain[i].nextPos = PosArray[i];
		chain[i + 1].prevPos = -PosArray[i];
	}
	chain[0].prevPos = 0;
	chain[n - 1].type = typeArray[n - 1];
	chain[n - 1].nextPos = 0;
	chain[n - 1].prevPos = -PosArray[n - 2];
	return 1;
}
Node* copychain(Node* chain) {
	Node* tmp = (Node*)malloc(LENGTH*sizeof(Node));
	if (tmp == NULL) {
		cout << "copychain malloc fail";
		exit(1);
	}
	for (int i = 0; i < LENGTH; i++) {
		tmp[i].type = chain[i].type;
		tmp[i].nextPos = chain[i].nextPos;
		tmp[i].prevPos = chain[i].prevPos;
	}
	return tmp;
}
char transPos(int i) {//翻译一下Pos
	switch (i) {
	case 1:return 'U';
	case -1:return 'D';
	case LENGTH:return 'L';
	case -LENGTH:return 'R';
	case LENGTH*LENGTH:return 'F';
	case -LENGTH*LENGTH:return 'B';
	case 0:return 'T';
	default:
		return 'E';
		break;
	}
}
void print(Node* chain) {//打印链
	for (int i = 0; i < LENGTH; i++) {
		cout << "(" << chain[i].type << "," << transPos(chain[i].nextPos) << ")" << "\t";
	}
	cout << endl;
}
void printInv(Node* chain) {//逆序打印链
	for (int i = LENGTH - 1; i >= 0; i--) {
		cout << "(" << chain[i].type << "," << transPos(chain[i].prevPos) << ")" << "\t";
	}
	cout << endl;
}
int randPos() {//产生随机Pos
	int r = rand() % 6;
	switch (r) {
	case 0:r = UP; break;
	case 1:r = DOWN; break;
	case 2:r = LEFT; break;
	case 3:r = RIGHT; break;
	case 4:r = FRONT; break;
	case 5:r = BEHIND; break;
	}
	return r;
}
/*int collision(Node* chain, int i) {//已被整合到Energy函数
switch (i)
{
case 0:{
int sum = 0;
for (int j = i; j < LENGTH; j++) {
sum += chain[j].nextPos;
if (sum == 0)
return 1;
}
return 0;
}
case LENGTH-1: {
int sum = 0;
for (int j = i; j >= 0; j--) {
sum += chain[j].prevPos;
if (sum == 0)
return 1;
}
return 0;
}
default: {
int sum = 0;
for (int j = i; j >= 0; j--) {
sum += chain[j].prevPos;
if (sum == 0)
return 1;
}
sum = 0;
for (int j = i; j < LENGTH; j++) {
sum += chain[j].nextPos;
if (sum == 0)
return 1;
}
return 0;
}
}
}*/
int Qube(Node* chain) {
	int xmax = -2*LENGTH*LENGTH, xmin = 2*LENGTH*LENGTH, ymax = -2*LENGTH*LENGTH, ymin = 2*LENGTH*LENGTH, zmax = -2*LENGTH*LENGTH, zmin = 2*LENGTH*LENGTH;
	int sumx = 0, sumy = 0, sumz = 0;
	for (int i = 0; i < LENGTH - 1; i++) {
		int nextPos = chain[i].nextPos;
		if (nextPos == UP || nextPos == DOWN) {
			sumx += nextPos;
			if (sumx < xmin)
				xmin = sumx;
			if (sumx > xmax)
				xmax = sumx;
		}
		if(nextPos == LEFT || nextPos == RIGHT) {
			sumy += nextPos;
			if (sumy < ymin)
				ymin = sumy;
			if (sumy > ymax)
				ymax = sumy;
		}
		if (nextPos == FRONT || nextPos == BEHIND) {
			sumz += nextPos;
			if (sumz < zmin)
				zmin = sumz;
			if (sumz > zmax)
				zmax = sumz;
		}
	}
	if (xmax - xmin == 2 && ymax - ymin == 2 && zmax - zmin == 2)
		return 1;
	else
		return 0;
}
int checkchain(Node* chain) {//检查整条链nextPos和prevPos之和是否为0
	int sum = 0;
	for (int i = 0; i < LENGTH; i++) {
		sum += chain[i].nextPos;
		sum += chain[i].prevPos;
	}
	if (sum != 0)
		return 1;
	else
		return 0;
}
int EnergyBetw(int a, int b) {//返回两个种类氨基酸的Energy
	if (a == b)
		return -3;
	else
		return -1;
}
int Energy(const Node chain[], int i) {//返回一个节点的Energy
	switch (i) {
	case 0:
	case 1: {
		int sum = 0, e = 0;
		for (int j = i; j < LENGTH; j++) {
			if ((sum == 1 || sum == -1 || sum == 27 || sum == -27 || sum == 729 || sum == -729) && abs(j - i)>1) {
				e += EnergyBetw(chain[i].type, chain[j].type);
			}
			sum += chain[j].nextPos;
			if (sum == 0) {
				e += 10000;//发生碰撞，能量很高
			}
		}
		return e;
		break;
	}
	case LENGTH - 1:
	case LENGTH - 2: {
		int sum = 0, e = 0;
		for (int j = i; j >= 0; j--) {
			if ((sum == 1 || sum == -1 || sum == 27 || sum == -27 || sum == 729 || sum == -729) && abs(j - i)>1) {
				e += EnergyBetw(chain[i].type, chain[j].type);
			}
			sum += chain[j].prevPos;
			if (sum == 0) {
				e += 10000;//发生碰撞，能量很高
			}
		}
		return e;
		break;
	}
	default: {
		int sum = 0, e = 0;
		for (int j = i; j >= 0; j--) {
			if ((sum == 1 || sum == -1 || sum == 27 || sum == -27 || sum == 729 || sum == -729) && abs(j - i)>1) {
				e += EnergyBetw(chain[i].type, chain[j].type);
			}
			sum += chain[j].prevPos;
			if (sum == 0) {
				e += 10000;//发生碰撞，能量很高
			}
		}
		sum = 0;
		for (int j = i; j < LENGTH; j++) {
			if ((sum == 1 || sum == -1 || sum == 27 || sum == -27 || sum == 729 || sum == -729) && abs(j - i)>1) {
				e += EnergyBetw(chain[i].type, chain[j].type);
			}
			sum += chain[j].nextPos;
			if (sum == 0) {
				e += 10000;//发生碰撞，能量很高
			}
		}
		return e;
	}
	}
}
int Energy(const Node chain[]) {//返回整条链的Energy，还可以优化减少一半时间
	int Esum = 0;
	for (int i = 0; i < LENGTH; i++) {
		Esum += Energy(chain, i);
	}
	return Esum / 2;
}
int CM(Node chain[], int i, double T) {//Corner Move
	int IEnergy = Energy(chain, i);
	int tmp;
	tmp = chain[i].nextPos;
	chain[i].nextPos = -chain[i].prevPos;
	chain[i].prevPos = -tmp;
	chain[i - 1].nextPos = -chain[i].prevPos;
	chain[i + 1].prevPos = -chain[i].nextPos;
	int dEnergy = Energy(chain, i) - IEnergy;
	if (dEnergy > 1000) {
		tmp = chain[i].nextPos;
		chain[i].nextPos = -chain[i].prevPos;
		chain[i].prevPos = -tmp;
		chain[i - 1].nextPos = -chain[i].prevPos;
		chain[i + 1].prevPos = -chain[i].nextPos;//发生碰撞，换回来
		/*if (checkchain(chain)) {
			cout << "CM Error" << endl;
			print(chain);
			exit(1);
		}*/
		return 0;
	}
	if (dEnergy < 0) {
		/*if (checkchain(chain)) {
			cout << "CM Error" << endl;
			print(chain);
			exit(1);
		}*/
		return dEnergy;
	}
	else {
		if (rand() / (double)RAND_MAX < exp(-dEnergy / T)) {
		/*	if (checkchain(chain)) {
				cout << "CM Error" << endl;
				print(chain);
				exit(1);
			}*/
			return dEnergy;
		}
		else {
			tmp = chain[i].nextPos;
			chain[i].nextPos = -chain[i].prevPos;
			chain[i].prevPos = -tmp;
			chain[i - 1].nextPos = -chain[i].prevPos;
			chain[i + 1].prevPos = -chain[i].nextPos;//大于概率，换回来	
			/*if (checkchain(chain)) {
				cout << "CM Error" << endl;
				print(chain);
				exit(1);
			}*/
			return 0;
		}
	}
}
int CSM(Node chain[], int i, double T) {//CrankShaft Move
	int direc = 757 - abs(chain[i].nextPos) - abs(chain[i].prevPos);//减法的结果就是剩下的那种Pos
	int tmp = chain[i].prevPos;
	int IEnergy = Energy(chain, i) + Energy(chain, i + 1);
	int r = rand() % 2;
	if (r == 0)
		r = -1;//r为0或-1的随机数
	chain[i].prevPos = r*direc;
	chain[i - 1].nextPos = -r*direc;
	chain[i + 1].nextPos = r*direc;
	chain[i + 2].prevPos = -r*direc;
	int dEnergy = Energy(chain, i) + Energy(chain, i + 1) - IEnergy;
	if (dEnergy > 1000) {
		chain[i].prevPos = tmp;
		chain[i - 1].nextPos = -tmp;
		chain[i + 1].nextPos = tmp;
		chain[i + 2].prevPos = -tmp;//发生碰撞，转回来
		/*if (checkchain(chain)) {
			cout << "CSM Error" << endl;
			print(chain);
			exit(1);
		}*/
		return 0;
	}
	if (dEnergy < 0) {
		/*if (checkchain(chain)) {
			cout << "CSM Error" << endl;
			print(chain);
			exit(1);
		}*/
		return dEnergy;
	}
	else {
		if (rand() / (double)RAND_MAX < exp(-dEnergy / T)) {
			/*if (checkchain(chain)) {
				cout << "CSM Error" << endl;
				print(chain);
				exit(1);
			}*/
			return dEnergy;
		}
		else {
			chain[i].prevPos = tmp;
			chain[i - 1].nextPos = -tmp;
			chain[i + 1].nextPos = tmp;
			chain[i + 2].prevPos = -tmp;//大于概率，转回来
			/*if (checkchain(chain)) {
				cout << "CSM Error" << endl;
				print(chain);
				exit(1);
			}*/
			return 0;
		}
	}
}
int EM(Node chain[], int i, double T) {//End Move
	int r;
	int IEnergy = Energy(chain, i);
	if (i == 0) {
		int tmp = chain[i].nextPos;
		do {
			r = randPos();
		} while (r == chain[i].nextPos || r == -chain[i].nextPos);
		chain[i].nextPos = r;
		chain[i + 1].prevPos = -r;
		int dEnergy = Energy(chain, i) - IEnergy;
		if (dEnergy > 1000) {
			chain[i].nextPos = tmp;
			chain[i + 1].prevPos = -tmp;//发生碰撞，转回来
			/*if (checkchain(chain)) {
				cout << "EM Error" << endl;
				print(chain);
				exit(1);
			}*/
			return 0;
		}
		if (dEnergy < 0) {
			/*if (checkchain(chain)) {
				cout << "EM Error" << endl;
				print(chain);
				exit(1);
			}*/
			return dEnergy;
		}
		else {
			if (rand() / (double)RAND_MAX < exp(-dEnergy / T)) {
				/*if (checkchain(chain)) {
					cout << "EM Error" << endl;
					print(chain);
					exit(1);
				}*/
				return dEnergy;
			}
			else {
				chain[i].nextPos = tmp;
				chain[i + 1].prevPos = -tmp;//大于概率，转回来
				/*if (checkchain(chain)) {
					cout << "EM Error" << endl;
					exit(1);
				}*/
				return 0;
			}
		}
	}
	if (i == LENGTH - 1) {
		int tmp = chain[i].prevPos;
		do {
			r = randPos();
		} while (r == chain[i].prevPos || r == -chain[i].prevPos);
		chain[i].prevPos = r;
		chain[i - 1].nextPos = -r;
		int dEnergy = Energy(chain, i) - IEnergy;
		if (dEnergy > 1000) {
			chain[i].prevPos = tmp;
			chain[i - 1].nextPos = -tmp;//发生碰撞，转回来
			/*if (checkchain(chain)) {
				cout << "EM Error" << endl;
				print(chain);
				exit(1);
			}*/
			return 0;
		}
		if (dEnergy < 0) {
			/*if (checkchain(chain)) {
				cout << "EM Error" << endl;
				print(chain);
				exit(1);
			}*/
			return dEnergy;
		}
		else {
			if (rand() / (double)RAND_MAX < exp(-dEnergy / T)) {
				/*if (checkchain(chain)) {
					cout << "EM Error" << endl;
					print(chain);
					exit(1);
				}*/
				return dEnergy;
			}
			else {
				chain[i].prevPos = tmp;
				chain[i - 1].nextPos = -tmp;//大于概率，转回来
				/*if (checkchain(chain)) {
					cout << "EM Error" << endl;
					print(chain);
					exit(1);
				}*/
				return 0;
			}
		}
	}
	print(chain);
	exit(3);//i不是1也不是LENGTH-1，出错
}
int judge(const Node chain[], int i) {//判断节点适用于哪种Move,0是EM，-1是CM，返回正数表示CSM的较小序号i，-2表示不动(在一条直线中间)
	if (i == 0 || i == LENGTH - 1) {
		return 0;
	}
	if (chain[i].nextPos == -chain[i].prevPos) {
		return -2;
	}
	else {
		int sum1 = chain[i].prevPos + chain[i + 1].prevPos + chain[i + 2].prevPos;
		int sum2 = chain[i].nextPos + chain[i - 1].nextPos + chain[i - 2].nextPos;
		if (sum1 == UP || sum1 == DOWN || sum1 == LEFT || sum1 == RIGHT || sum1 == FRONT || sum1 == BEHIND)
			return i;
		if (sum2 == UP || sum2 == DOWN || sum2 == LEFT || sum2 == RIGHT || sum2 == FRONT || sum2 == BEHIND)
			return i - 1;
		return  -1;
	}
}
int MCstep(Node chain[], double T) {//蒙特卡罗一步，未执行返回0;执行返回dEnergy;
	int i = rand() % LENGTH;
	int re = judge(chain, i);
	if (re > 0) {
		return CSM(chain, re, T);
	}
	else {
		switch (re) {
		case 0:return EM(chain, i, T); break;
		case -1:return CM(chain, i, T); break;
		case -2:return 0; break;
		default:exit(2);
		}
	}
}

double MeanE(long long int a[], int n) {//求能量分布期望
	double mean = 0;
	long long int num = 0;
	for (int i = 0; i < n; i++) {
		num += a[i];
		mean += -i*a[i];
	}
	if (num == 0)
		return 0;
	return mean / num;
}
double VarianceE(long long int a[], int n) {//求能量分布方差
	double mean = MeanE(a, n), sum = 0, num = 0;
	for (int i = 0; i < n; i++) {
		num += a[i];
		sum += a[i] * i*i;
	}
	return sum / num - mean*mean;
}


int main()
{
	srand(time(NULL));
	int typeArrayBest[LENGTH] = { 0,1,0,1,1,1,1,1,0,1,1,0,1,0,1,0,0,0,1,1,0,0,0,0,0,0,1 };//ABABBBBBABBABABAAABBAAAAAAB
	int typeArray[LENGTH] = { 0 };
	int PosArray[LENGTH - 1] = { DOWN,DOWN,RIGHT,RIGHT,UP,UP,UP,UP,LEFT,LEFT,LEFT,LEFT,DOWN,DOWN,DOWN,DOWN,DOWN,DOWN,RIGHT,RIGHT,RIGHT,RIGHT,RIGHT,RIGHT,UP,UP };
	int PosArrayBest[LENGTH - 1] = { RIGHT,FRONT,RIGHT,BEHIND,UP,UP,LEFT,FRONT,RIGHT,DOWN,LEFT,BEHIND,LEFT,UP,FRONT,FRONT,RIGHT,RIGHT,DOWN,LEFT,LEFT,BEHIND,DOWN,FRONT,RIGHT,RIGHT };

	
	Node chain[LENGTH];
	InitNode(chain, LENGTH, typeArrayBest);	
	Node chainBest[LENGTH];
	InitNode(chainBest, LENGTH, PosArray, typeArrayBest);

	/*
	printInv(chain);
	int Esum = 0;
	cout << Energy(chain, 5) << endl;
	cout << endl;
	printInv(chainBest);
	for (int i = 0; i < LENGTH; i++) {
		cout << Energy(chainBest, i) << " ";
		Esum += Energy(chainBest, i);
	}
	cout << endl << Esum / 2 << endl;
	cout << judge(chainBest, 1) << endl;

	cout << endl;
	*/
	//以上为测试部分

	
	//Fold time test
	/*int currentEnergy = Energy(chainBest);
	cout << "InitEnergy=\t" << currentEnergy << endl;

	ofstream f1("D:\\Output\\Output T=2.8-3.0 foldtime.txt");
	int sum = 0;
	int finish = 0, NumUnfinish = 0;
	long long int i = 0;
	long long int N[NUMTEST] = { 0 };

	cout << endl << "Begin" << endl;
	double T;
	for (T = 2.8; T <= 3.01; T += 0.2) {
		NumUnfinish = 0;
		for (int j = 0; j < NUMTEST; j++) {
			InitNode(chainBest, LENGTH, PosArray, typeArrayBest);
			currentEnergy = Energy(chaintest);
			finish = 0;
			i = 0;
			N[j] = 0;
			while (i < 1000000000) {
				if (currentEnergy == -84) {
					N[j] = i;
					finish = 1;
					break;
				}
				currentEnergy += MCstep(chainBest, T);
				i++;
			}
			if (finish == 0) {
				NumUnfinish++;
				N[j] = -1;
			}
		}
		f1 << "T = " << T << endl;
		for (int i = 0; i < NUMTEST; i++) {
			f1 << N[i] << endl;
		}
		f1 << "T=" << T << "\tNumUnfinish=" << NumUnfinish << endl << endl;
		cout << "T=" << T << "\tfinished" << endl;
	}	*/
	
	
	//Stablility test
	
	/*Node chainFolded[LENGTH];
	InitNode(chainFolded, LENGTH, PosArrayBest, typeArrayBest);

	ofstream f2("D:\\Output\\Output folded test.txt");
	long long int E[85] = { 0 };
	int currentEnergy = Energy(chainFolded);
	cout << "InitEnergy=\t" << currentEnergy << endl;

	double T;
	cout << endl << "Begin" << endl;
	for (T = 0.8; T <= 3.01; T += 0.2) {
		InitNode(chainFolded, LENGTH, PosArrayBest, typeArrayBest);
		currentEnergy = Energy(chaintest);
		for (int i = 0; i < 85; i++) {
			E[i] = 0;
		}
		for (long long int i = 0; i < 10000000; i++) {
			E[-currentEnergy]++;
			currentEnergy += MCstep(chainFolded, T);
		}
		f2 << "T = " << T << endl;
		for (int i = 0; i < 85; i++) {
			f2 << E[i] << endl;
		}
		f2 << endl;
	}*/

	//Energy histogram

	/*Node chaintest[LENGTH];
	InitNode(chaintest, LENGTH, PosArray, typeArrayBest);

	ofstream f2("D:\\Output\\Output Energy Histograms T=2-5.txt");
	
	double T[4] = { 2.00,2.51,3.15,5.00 };
	int currentEnergy = Energy(chaintest);
	cout << "InitEnergy=\t" << currentEnergy << endl;

	cout << endl << "Begin" << endl;
	for (int j = 0; j < 4; j++) {
		InitNode(chaintest, LENGTH, PosArray, typeArrayBest);
		currentEnergy = Energy(chaintest);
		long long int E[85] = { 0 };
		for (int i = 0; i < 85; i++) {
			E[i] = 0;
		}
		for (long long int i = 0; i < 1000000000; i++) {
			E[-currentEnergy]++;
			currentEnergy += MCstep(chaintest, T[j]);
		}
		f2 << "T = " << T[j] << endl;
		for (int i = 0; i < 85; i++) {
			f2 << E[i] << endl;
		}
		f2 << endl;
		cout << "T = " << T[j] << "\tfinished" << endl;
	}*/


	//Average Energy of the first τ steps
	Node chainME[LENGTH];
	int currentEnergytest = 0;
	double meanEArr[tINACF/1000] = { 0 };
	for (int i = 0; i < tINACF /1000; i++) {
		meanEArr[i] = 0;
	}
	long long int Etest[85] = { 0 };
	for (int i = 0; i < 85; i++) {
		Etest[i] = 0;
	}

	for (int j = 0; j < 100; j++) {
		InitNode(chainME, LENGTH, PosArray, typeArrayBest);
		currentEnergytest = Energy(chainME);
		for (int i = 0; i < 85; i++) {
			Etest[i] = 0;
		}
		for (int i = 0; i < tINACF; i++) {
			if (i % 1000 == 0) {
				meanEArr[i / 1000] += MeanE(Etest, 85);
			}
			Etest[-currentEnergytest]++;
			currentEnergytest += MCstep(chainME, 1.58);
		}
	}
	for (int i = 0; i < tINACF / 1000; i++) {
		meanEArr[i] = meanEArr[i] / 100;
	}

	//autocorrection function
	Node chaintest2[LENGTH];
	InitNode(chaintest2, LENGTH, PosArray, typeArrayBest);
	ofstream f2("D:\\Output\\Output Autocorrection Function T=1.58 times=1000.txt");	
	double T[4] = { 1.58,1.30,1.45,1.80 };
	
	int currentEnergy = Energy(chaintest2);
	cout << "InitEnergy=\t" << currentEnergy << endl;

	cout << endl << "Begin" << endl;
	for (int j = 0; j < 1; j++) {
		double ACF[tINACF / 1000] = { 0 };
		for (int i = 0; i <tINACF / 1000; i++)
			ACF[i] = 0.0;
		for (int k = 0; k <ACFSTEP; k++) {
			InitNode(chaintest2, LENGTH, PosArray, typeArrayBest);
			currentEnergy = Energy(chaintest2);
			long long int Etest2[85] = { 0 };
			for (int i = 0; i < 85; i++) {
				Etest2[i] = 0;
			}
			for (long long int i = 0; i < τ; i++) {
				Etest2[-currentEnergy]++;
				currentEnergy += MCstep(chaintest2, T[j]);
			}//先runτ次
			double initE = Energy(chaintest2);
			double C0 = 0;
			for (int i = 0; i < 85; i++) {
				C0 += Etest2[i] *i*i/ τ;
			}
			cout << initE << "\t,\t" << C0 << endl;
			long long int E[85] = { 0 };
			for (int i = 0; i < 85; i++) {
				E[i] = 0;
			}
			for (long long int i = 0; i < tINACF;i++) {
				E[-currentEnergy]++;
				currentEnergy += MCstep(chaintest2, T[j]);
				
				if (i % 1000 == 0) {
					ACF[i / 1000] += ((initE * currentEnergy - meanEArr[i / 1000] * meanEArr[i / 1000]) / C0);
				}
			}
		}
		f2 << "T = " << T[j] << endl;
		for (int i = 0; i < tINACF/1000; i++)
			ACF[i] = ACF[i] / ACFSTEP;
		for (int i = 0; i < tINACF/1000; i++)
			f2 << ACF[i] << endl;
		f2 << endl;
		cout << "T = " << T[j] << "\tfinished" << endl;
	}

	return 0;
}

