#pragma once
#define LENGTH 27
#define UP 1
#define DOWN -1
#define LEFT LENGTH
#define RIGHT -LENGTH
#define FRONT LENGTH*LENGTH
#define BEHIND -LENGTH*LENGTH
using namespace std;

struct Node {
	int type;//0��ʾA���࣬1��ʾB����
	int nextPos;//���������27����1,-1 27,-27 729,-729��Ӧ��,�� ��,�� ǰ,��0��ʾ��
	int prevPos;
}*Chain;

int InitNode(Node chain[], int n, int typeArray[]);
int InitNode(Node chain[], int n, int PosArray[], int typeArray[]);
Node* copychain(Node* chain);
char transPos(int i);
void print(Node* chain);
void printInv(Node* chain);
int randPos();

int Qube(Node* chain);
int checkchain(Node* chain);
int EnergyBetw(int a, int b);
int Energy(const Node chain[], int i);
int Energy(const Node chain[]);
int CM(Node chain[], int i, double T);
int CSM(Node chain[], int i, double T);
int EM(Node chain[], int i, double T);
int judge(const Node chain[], int i);
int MCstep(Node chain[], double T);

double Mean(long long int a[], int n);
double StandardDeviation(long long int a[], int n);



