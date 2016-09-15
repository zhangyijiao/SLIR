#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mt19937ar.c"

struct Link
{
	int num;
	struct Link *next;
};

void AddNode(struct Link *net, int i) 
{
	struct Link *p;
	p = (struct Link *)malloc(sizeof(struct Link));
	p->next = net->next;
	p->num = i;
	net->next = p;
}

void AddLink(struct Link *net, int i, int j)
{
	AddNode(&net[i], j);
	AddNode(&net[j], i);
}

void CreateNetwork(struct Link *net, int n, double probability)
{
	int i, j;
	double ran;
	for(i = 0; i < n; i++) {
		for(j = i + 1; j < n; j++) {
			ran = genrand_real3();
			if(ran < probability) {
				AddLink(net, i, j);
			}
		}
	}
}

//state: 'S'--'0', 'A'--'10', 'B'--'1', 'I'--'11', 'R'--'2'
int Spread(struct Link *net, int n, double a1, double a2) 
{
	int i, number_A, number_B, number_I, number_R, ran, temp, step;
	struct Link *p;
	double c1, c2, k1, k2, k3;
	int *degree = (int *)malloc(n*sizeof(int));
	int *state = (int *)malloc(n*sizeof(int));
	int *temp_state = (int *)malloc(n*sizeof(int));
	//calculate nodes' degree and choose a seed whose degree isn't zero
	for(i = 0; i < n; i++) {
		temp = 0;
		p = &net[i];
		while(p->next != NULL) {
			p = p->next;
			temp += 1;
		}
		degree[i] = temp;
	}
	ran = (int)(n*genrand_real1());
	while(degree[ran] == 0) {
		ran = (int)(n*genrand_real1());
	}
	//initial nodes' state and start to spread
	for(i = 0; i < n; i++) {
		if(i == ran) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_A = 0;
	number_B = 0;
	number_R = 0;
	step = 0;
	while((number_I > 0) || ((number_A > 0) && (number_B > 0))) {
		for(i = 0; i < n; i++) {
			//node is S
			if(state[i] == 0) {
				p = &net[i];
				k1 = 0.;
				k2 = 0.;
				k3 = 0.;
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 10) k1 += 1;
					else if(state[p->num] == 1) k2 += 1;
					else if(state[p->num] == 11) k3 +=1;
				}
				c1 = 1 - pow((1 - a1), k1)*pow((1 - a2), k3);
				c2 = 1 - pow((1 - a1), k2)*pow((1 - a2), k3);
				if(genrand_real2() < c1) {
					temp_state[i] = 10;
					number_A += 1;
					if(genrand_real2() < c2) {
						temp_state[i] = 11;
						number_A -= 1;
						number_I += 1;
					}
				}
				else {
					if(genrand_real2() < c2) {
						temp_state[i] = 1;
						number_B += 1;
					}
					else temp_state[i] = 0;
				}
			}
			//node is A
			else if(state[i] == 10) {
				p = &net[i];
				k2 = 0.;
				k3 = 0.;
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 1) k2 +=1;
					else if(state[p->num] == 11) k3 +=1;
				}
				c2 = 1 - pow((1 - a1), k2)*pow((1 - a2), k3);
				if(genrand_real2() < c2) {
					temp_state[i] = 11;
					number_A -= 1;
					number_I += 1;
				}
				else temp_state[i] = 10;
			}
			//node is B
			else if(state[i] == 1) {
				p = &net[i];
				k1 = 0.;
				k3 = 0.;
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 10) k1 +=1;
					else if(state[p->num] == 11) k3 +=1;
				}
				c1 = 1 - pow((1 - a1), k1)*pow((1 - a2), k3);
				if(genrand_real2() < c1) {
					temp_state[i] = 11;
					number_B -= 1;
					number_I += 1;
				}
				else temp_state[i] = 1;
			}
			//node is I
			else if(state[i] == 11) {
				temp_state[i] = 2;
				number_I -= 1;
				number_R += 1;
			}
			//node is R
			else {
				if(state[i] != 2) printf("error~~~~~~~\n");
				temp_state[i] = 2;
			}
		}
		for(i = 0; i < n; i++) {
			state[i] = temp_state[i];
		}
		step += 1;
		if(step > 100) break;
	}
	free(degree);
	free(state);
	free(temp_state);
	return number_R;
}

int main(int argc, char *argv[])
{
	int i, j, n = 10000, number_R, temp, total_degree = 0;
	double a1, a2;
	double probability = 4./(n -1);
	struct Link *net = (struct Link*)malloc(n*sizeof(struct Link));
	struct Link *p, *pp;
	char file_name[25];
	int b1, b2;
	if (argc > 2) {
		a1 = atof(argv[1]);
		a2 = atof(argv[2]);
	} else {
		printf("no parameter input\n");
		exit(1);
	}

	b1 = (int)(a1*100);
	b2 = (int)(a2*100);
	sprintf(file_name, "SLIR_random_%d_%d.txt", b1, b2);
	FILE *fp = fopen(file_name, "w");
	fprintf(fp, "#a1=%.2f\n", a1);
	fprintf(fp, "#a2=%.2f\n", a2);

	unsigned long idum = (unsigned long)time(NULL);
	init_genrand(idum);
	for(i = 0; i < n; i++) {
		net[i].num = i;
		net[i].next = NULL;
	}
	for(i = 0; i < 100; i++) {
		CreateNetwork(net, n, probability);
		for(j = 0; j < 100; j++) {
			number_R = Spread(net, n, a1, a2);
			fprintf(fp, "%d\n", number_R);
		}
		
		for(j = 0; j < n; j++) {
			p = &net[j];
			while(p->next != NULL) {
				pp = p->next;
				p->next = pp->next;
				free(pp);
			}
		}
	}
	fclose(fp);
	printf("meow\n");
	return 0;
}