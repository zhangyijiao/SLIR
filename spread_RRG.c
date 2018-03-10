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

void FreeNetwork(struct Link *net, int n)
{
	int i;
	struct Link *p, *pp;
	for(i = 0; i < n; i++) {
		p = &net[i];
		while(p->next != NULL) {
			pp = p->next;
			p->next = pp->next;
			free(pp);
		}
	}
}

int Suitable(struct Link *net, int i, int j)
{
    struct Link *p;
    int index = 0;
    p = &net[i];
    if(i == j) {
        index = 1;
    }
    while(p->next != NULL) {
        p = p->next;
        if(p->num == j){
            index = 1;
            break;
        }
    }
    return index;
}

void RegularRandomG(struct Link *net, int n, int k)
{
    int i, j, temp, node, number_node, temp_number_node;
    int *node_list, *temp_node_list;
    temp_number_node = 1;
    while(temp_number_node != 0) {
        FreeNetwork(net, n);
        number_node = n*k;
        //generate a node list according to the nodes degree
        node_list = (int *)malloc(number_node*sizeof(int));
        temp_node_list = (int *)malloc(number_node*sizeof(int));
        temp = 0;
        for(i = 0; i < n; i++) {
            for(j = 0; j < k; j++) {
                node_list[temp] = i;
                temp++;
            }
        }    
        //connect the suitable links
        while(number_node > 0) {
            //shuffle the list
            for(i = number_node - 1; i > 0; i--){
                node = (int)(i*genrand_real3());
                temp = node_list[i];
                node_list[i] = node_list[node];
                node_list[node] = temp;
            }
            temp_number_node = 0;
            for(i = 0; i < number_node; i+=2) {
                if(Suitable(net, node_list[i], node_list[i+1]) == 0){
                    AddLink(net, node_list[i], node_list[i+1]);
                    // number_node -= 2;
                }
                else {
                    temp_node_list[temp_number_node] = node_list[i];
                    temp_node_list[temp_number_node+1] = node_list[i+1];
                    temp_number_node += 2;
                }
            }
            if(number_node == temp_number_node) break;
            number_node = temp_number_node;
            free(node_list);
            node_list = (int *)malloc(number_node*sizeof(int));
            for(i = 0; i < number_node; i++) {
                node_list[i] = temp_node_list[i];
            }
            free(temp_node_list);
            temp_node_list = (int *)malloc(number_node*sizeof(int));
        }
    }    
}

//state: 'S'--'0', 'L'--'1', 'I'--'11', 'R'--'2'
//for regular graph such as lattice and RRG 
int SLIRStaticLattice(struct Link *net, int n, double a)
{
	int i, number_I, number_RL, number_R, ran_node, temp, max_degree;
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate;
	struct Link *p;
	double ran;
	max_degree = 6;
	ran_node = (int)(n*genrand_real3());
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	number_RL = 0;
	for(i = 0; i < n; i++) {
		if(state[i] == 11) rate[i] = 1;
		else {
			p = &net[i];
			temp = 0;
			while(p->next != NULL) {
				p = p->next;
				if(state[p->num] == 11) {
					temp += 1;
				}
			}
			rate[i] = temp*a*(2 - a);
		}
	}
	if(a*max_degree*(2 - a) > 1.0) max_rate = max_degree*a*(2 - a);
	else max_rate = 1.0;
	while(number_I > 0) {
		ran_node = (int)(n*genrand_real3());
		while(rate[ran_node] < 0.00001) {
			ran_node = (int)(n*genrand_real3());
		}
		if(genrand_real3() < rate[ran_node]/max_rate) {
			if(state[ran_node] == 11) {
				state[ran_node] = 2;
				rate[ran_node] = 0;
				number_I -= 1;
				number_R += 1;
				p = &net[ran_node];
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 0) rate[p->num] -= a*(2 - a);
					else if(state[p->num] == 1) rate[p->num] -= a;
				}
			}
			else if(state[ran_node] == 0) {
				ran = genrand_real3();
				if(ran < 2*(1 - a)/(2 - a)) {
					state[ran_node] = 1;
					p = &net[ran_node];
					temp = 0;
					while(p->next != NULL) {
						p = p->next;
						if(state[p->num] == 11) temp += 1;
					}
					rate[ran_node] = temp*a;
				}
				else {
					state[ran_node] = 11;
					rate[ran_node] = 1;
					number_I += 1;
					p = &net[ran_node];
					while(p->next != NULL) {
						p = p->next;
						if(state[p->num] == 0) rate[p->num] += a*(2 - a);
						else if(state[p->num] == 1) rate[p->num] += a;
					}
				}
			}
			else if(state[ran_node] == 1) {
				state[ran_node] = 11;
				rate[ran_node] = 1;
				number_I += 1;
				number_RL += 1;
				p = &net[ran_node];
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 0) rate[p->num] += a*(2 - a);
					else if(state[p->num] == 1) rate[p->num] += a;
				}
			}
		}
	}
	// printf("number_R:%d\n", number_R);
	free(state);
	free(rate);
	return number_RL;
}

int SLIRAnnealLattice(int n, double a, int k)
{
	int i, j, number_I, number_L, number_R, ran_node;
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate;
	struct Link *p;
	int ran, k1;
	double ran_num;
	int *neighbor = (int *)malloc(k*sizeof(int));
	int *label = (int *)malloc(n*sizeof(int));
	for(i = 0; i < n; i++) {
		label[i] = 0;
	}
	ran_node = (int)(n*genrand_real3());
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		rate[i] = 0;
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	number_L = 0;
	if(a*k*(2 - a) > 1.0) max_rate = k*a*(2 - a);
	else max_rate = 1.0;
	while(number_I > 0) {
		ran_node = (int)(n*genrand_real3());
		while(state[ran_node] == 2) {
			ran_node = (int)(n*genrand_real3());
		}
		if(state[ran_node] == 11) rate[ran_node] = 1;
		else if((state[ran_node] == 0) || (state[ran_node] == 1)) {
			k1 = 0;
			//choose neighbor
			for(j = 0; j < k; j++) {
				ran = (int)(n*genrand_real3());
				while((label[ran] == 1) || (ran == ran_node)) {
					ran = (int)(n*genrand_real3());
				}
				neighbor[j] = ran;
				label[ran] = 1;
				if(state[ran] == 11) k1 += 1;
			}
			if(state[ran_node] == 0) rate[ran_node] = k1*a*(2-a);
			else rate[ran_node] = k1*a;
		}
		if(genrand_real3() < rate[ran_node]/max_rate) {
			if(state[ran_node] == 11) {
				state[ran_node] = 2;
				number_I -= 1;
				number_R += 1;
			}
			else if(state[ran_node] == 0) {
				ran_num = genrand_real3();
				if(ran_num < 2*(1 - a)/(2 - a)) {
					state[ran_node] = 1;
					number_L += 1;
				}
				else {
					state[ran_node] = 11;
					number_I += 1;
				}
			}
			else if(state[ran_node] == 1) {
				state[ran_node] = 11;
				number_I += 1;
				number_L -= 1;
			}
		}
		rate[ran_node] = 0;
		for(j = 0; j < k; j++) {
			label[neighbor[j]] = 0;
		}
	}
	// printf("number_R:%d\n", number_R);
	free(state);
	free(rate);
	free(neighbor);
	free(label);
	return number_R;
}

int SLIRStatictimeOnce(struct Link *net, int n, double a)
{
	int i, number_I, number_L, number_R, ran_node, temp, max_degree;
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate, total_rate;
	struct Link *p;
	double ran, t = 0.;
	int temp_t = 50;
	int b;
	b = (int)(a*100);
	char file_name1[20];
	char file_name2[20];
	char file_name3[20];
	sprintf(file_name1, "SLIR_staticR_%d.txt", b);
	sprintf(file_name2, "SLIR_staticL_%d.txt", b);
	sprintf(file_name3, "SLIR_staticI_%d.txt", b);
	FILE *fp1 = fopen(file_name1, "w");
	FILE *fp2 = fopen(file_name2, "w");
	FILE *fp3 = fopen(file_name3, "w");
	max_degree = 6;
	ran_node = (int)(n*genrand_real3());
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	number_L = 0;
	total_rate = 0;
	for(i = 0; i < n; i++) {
		if(state[i] == 11) rate[i] = 1;
		else {
			p = &net[i];
			temp = 0;
			while(p->next != NULL) {
				p = p->next;
				if(state[p->num] == 11) {
					temp += 1;
				}
			}
			rate[i] = temp*a*(2 - a);
		}
		total_rate += rate[i];
	}
	if(a*max_degree*(2 - a) > 1.0) max_rate = max_degree*a*(2 - a);
	else max_rate = 1.0;
	while(number_I > 0) {
		ran_node = (int)(n*genrand_real3());
		while(rate[ran_node] < 0.00001) {
			ran_node = (int)(n*genrand_real3());
		}
		if(genrand_real3() < rate[ran_node]/max_rate) {
			//record result every 0.05s
			if((int)(t/0.05 + 0.000001) == temp_t) { 
				fprintf(fp1, "%d %d\n", temp_t, number_R);
				fprintf(fp2, "%d %d\n", temp_t, number_L);
				fprintf(fp3, "%d %d\n", temp_t, number_I);
				temp_t += 1;
			}
			t += 1.0/total_rate; 
			if(state[ran_node] == 11) {
				state[ran_node] = 2;
				total_rate -= rate[ran_node];
				rate[ran_node] = 0;
				number_I -= 1;
				number_R += 1;
				p = &net[ran_node];
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 0) {
						rate[p->num] -= a*(2 - a);
						total_rate -= a*(2 - a);
					}
					else if(state[p->num] == 1) {
						rate[p->num] -= a;
						total_rate -= a;
					} 
				}
			}
			else if(state[ran_node] == 0) {
				ran = genrand_real3();
				if(ran < 2*(1 - a)/(2 - a)) {
					state[ran_node] = 1;
					number_L += 1;
					p = &net[ran_node];
					temp = 0;
					while(p->next != NULL) {
						p = p->next;
						if(state[p->num] == 11) temp += 1;
					}
					rate[ran_node] = temp*a;
					total_rate += temp*a;
				}
				else {
					state[ran_node] = 11;
					rate[ran_node] = 1;
					total_rate += 1;
					number_I += 1;
					p = &net[ran_node];
					while(p->next != NULL) {
						p = p->next;
						if(state[p->num] == 0) {
							rate[p->num] += a*(2 - a);
							total_rate += a*(2 - a);
						} 
						else if(state[p->num] == 1) {
							rate[p->num] += a;
							total_rate += a;
						} 
					}
				}
			}
			else if(state[ran_node] == 1) {
				state[ran_node] = 11;
				rate[ran_node] = 1;
				total_rate += 1;
				number_I += 1;
				number_L -= 1;
				p = &net[ran_node];
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 0) {
						rate[p->num] += a*(2 - a);
						total_rate += a*(2 - a);
					} 
					else if(state[p->num] == 1) {
						rate[p->num] += a;
						total_rate += a;
					} 
				}
			}
		}
	}
	// printf("number_R:%d \n", number_R);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	free(state);
	free(rate);
	return number_R;
}

int SLIRAnnealtimeOnce(int n, double a, int k)
{
	int i, j, number_I, number_L, number_R, ran_node;
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate;
	int ran, k1;
	int *neighbor = (int *)malloc(k*sizeof(int));
	int *label = (int *)malloc(n*sizeof(int));
	int temp_t = 50;
	double total_rate, t = 0.;
	int b;
	b = (int)(a*100);
	char file_name1[20];
	char file_name2[20];
	char file_name3[20];
	sprintf(file_name1, "SLIR_annealR_%d.txt", b);
	sprintf(file_name2, "SLIR_annealL_%d.txt", b);
	sprintf(file_name3, "SLIR_annealI_%d.txt", b);
	FILE *fp1 = fopen(file_name1, "w");
	FILE *fp2 = fopen(file_name2, "w");
	FILE *fp3 = fopen(file_name3, "w");
	ran_node = (int)(n*genrand_real3());
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		label[i] = 0;
		rate[i] = 0;
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	number_L = 0;
	if(a*k*(2 - a) > 1.0) max_rate = k*a*(2 - a);
	else max_rate = 1.0;
	while(number_I > 0) {
		ran_node = (int)(n*genrand_real3());
		while(state[ran_node] == 2) {
			ran_node = (int)(n*genrand_real3());
		}
		if(state[ran_node] == 11) rate[ran_node] = 1;
		else if((state[ran_node] == 0) || (state[ran_node] == 1)) {
			k1 = 0;
			//choose neighbor
			for(j = 0; j < k; j++) {
				ran = (int)(n*genrand_real3());
				while((label[ran] == 1) || (ran == ran_node)) {
					ran = (int)(n*genrand_real3());
				}
				neighbor[j] = ran;
				label[ran] = 1;
				if(state[ran] == 11) k1 += 1;
			}
			if(state[ran_node] == 0) rate[ran_node] = k1*a*(2-a);
			else rate[ran_node] = k1*a;
		}
		if(genrand_real3() < rate[ran_node]/max_rate) {
			if((int)(t/0.1 + 0.000001) == temp_t) {
				fprintf(fp1, "%d %d\n", temp_t, number_R);
				fprintf(fp2, "%d %d\n", temp_t, number_L);
				fprintf(fp3, "%d %d\n", temp_t, number_I);
				// printf("temp_t:%d  time: %f  R:%d \n", temp_t, t, number_R);
				temp_t += 1;
			}
			total_rate = number_I + number_L*k*number_I*a/n + (n - number_I - number_L - number_R)*k*number_I*a*(2-a)/n;
			t += 1/total_rate;
			if(state[ran_node] == 11) {
				state[ran_node] = 2;
				number_I -= 1;
				number_R += 1;
			}
			else if(state[ran_node] == 0) {
				if(genrand_real3() < 2*(1 - a)/(2 - a)) {
					state[ran_node] = 1;
					number_L += 1;
				}
				else {
					state[ran_node] = 11;
					number_I += 1;
				}
			}
			else if(state[ran_node] == 1) {
				state[ran_node] = 11;
				number_I += 1;
				number_L -= 1;
			}
		}
		rate[ran_node] = 0;
		for(j = 0; j < k; j++) {
			label[neighbor[j]] = 0;
		}
	}
	// printf("number_R:%d time:%f\n", number_R, t);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	free(state);
	free(rate);
	free(neighbor);
	free(label);
	return number_R;
}

int SLIRStatictime(struct Link *net, int n, double a)
{
	int i, number_I, number_L, number_R, ran_node, temp, max_degree;
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate, total_rate;
	struct Link *p;
	double ran, t = 0.;
	int temp_t = 50;
	char file_name1[20];
	char file_name2[20];
	char file_name3[20];
	max_degree = 6;
	ran_node = (int)(n*genrand_real3());
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	number_L = 0;
	total_rate = 0;
	for(i = 0; i < n; i++) {
		if(state[i] == 11) rate[i] = 1;
		else {
			p = &net[i];
			temp = 0;
			while(p->next != NULL) {
				p = p->next;
				if(state[p->num] == 11) {
					temp += 1;
				}
			}
			rate[i] = temp*a*(2 - a);
		}
		total_rate += rate[i];
	}
	if(a*max_degree*(2 - a) > 1.0) max_rate = max_degree*a*(2 - a);
	else max_rate = 1.0;
	while(number_I > 0) {
		ran_node = (int)(n*genrand_real3());
		while(rate[ran_node] < 0.00001) {
			ran_node = (int)(n*genrand_real3());
		}
		// t += 1.0/total_rate; //73
		if(genrand_real3() < rate[ran_node]/max_rate) {
			//record result every 0.05s
			if((int)(t/0.05 + 0.000001) == temp_t) { 
				sprintf(file_name1,"SLIR_staticR_%d.txt", temp_t);
				FILE *fp1 = fopen(file_name1, "a");
				fprintf(fp1, "%d\n", number_R);
				fclose(fp1);
				sprintf(file_name2,"SLIR_staticL_%d.txt", temp_t);
				FILE *fp2 = fopen(file_name2, "a");
				fprintf(fp2, "%d\n", number_L);
				fclose(fp2);
				sprintf(file_name3,"SLIR_staticI_%d.txt", temp_t);
				FILE *fp3 = fopen(file_name3, "a");
				fprintf(fp3, "%d\n", number_I);
				fclose(fp3);
				temp_t += 1;
			}
			t += 1.0/total_rate; //12
			if(state[ran_node] == 11) {
				state[ran_node] = 2;
				total_rate -= rate[ran_node];
				rate[ran_node] = 0;
				number_I -= 1;
				number_R += 1;
				p = &net[ran_node];
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 0) {
						rate[p->num] -= a*(2 - a);
						total_rate -= a*(2 - a);
					}
					else if(state[p->num] == 1) {
						rate[p->num] -= a;
						total_rate -= a;
					} 
				}
			}
			else if(state[ran_node] == 0) {
				ran = genrand_real3();
				if(ran < 2*(1 - a)/(2 - a)) {
					state[ran_node] = 1;
					number_L += 1;
					p = &net[ran_node];
					temp = 0;
					while(p->next != NULL) {
						p = p->next;
						if(state[p->num] == 11) temp += 1;
					}
					rate[ran_node] = temp*a;
					total_rate += temp*a;
				}
				else {
					state[ran_node] = 11;
					rate[ran_node] = 1;
					total_rate += 1;
					number_I += 1;
					p = &net[ran_node];
					while(p->next != NULL) {
						p = p->next;
						if(state[p->num] == 0) {
							rate[p->num] += a*(2 - a);
							total_rate += a*(2 - a);
						} 
						else if(state[p->num] == 1) {
							rate[p->num] += a;
							total_rate += a;
						} 
					}
				}
			}
			else if(state[ran_node] == 1) {
				state[ran_node] = 11;
				rate[ran_node] = 1;
				total_rate += 1;
				number_I += 1;
				number_L -= 1;
				p = &net[ran_node];
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 0) {
						rate[p->num] += a*(2 - a);
						total_rate += a*(2 - a);
					} 
					else if(state[p->num] == 1) {
						rate[p->num] += a;
						total_rate += a;
					} 
				}
			}
		}
	}
	// printf("number_R:%d \n", number_R);
	free(state);
	free(rate);
	return number_R;
}

int SLIRAnnealtime(int n, double a, int k)
{
	int i, j, number_I, number_L, number_R, ran_node;
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate;
	struct Link *p;
	int ran, k1;
	double ran_num;
	int *neighbor = (int *)malloc(k*sizeof(int));
	int *label = (int *)malloc(n*sizeof(int));
	int temp_t = 50;
	double total_rate, t = 0.;
	char file_name1[20];
	char file_name2[20];
	char file_name3[20];
	// FILE *fp = fopen("SLIR_anneal_time.txt", "w");
	for(i = 0; i < n; i++) {
		label[i] = 0;
	}
	ran_node = (int)(n*genrand_real3());
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		rate[i] = 0;
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	number_L = 0;
	if(a*k*(2 - a) > 1.0) max_rate = k*a*(2 - a);
	else max_rate = 1.0;
	while(number_I > 0) {
		ran_node = (int)(n*genrand_real3());
		while(state[ran_node] == 2) {
			ran_node = (int)(n*genrand_real3());
		}
		if(state[ran_node] == 11) rate[ran_node] = 1;
		else if((state[ran_node] == 0) || (state[ran_node] == 1)) {
			k1 = 0;
			//choose neighbor
			for(j = 0; j < k; j++) {
				ran = (int)(n*genrand_real3());
				while((label[ran] == 1) || (ran == ran_node)) {
					ran = (int)(n*genrand_real3());
				}
				neighbor[j] = ran;
				label[ran] = 1;
				if(state[ran] == 11) k1 += 1;
			}
			if(state[ran_node] == 0) rate[ran_node] = k1*a*(2-a);
			else rate[ran_node] = k1*a;
		}
		if(genrand_real3() < rate[ran_node]/max_rate) {
			if((int)(t/0.1 + 0.000001) == temp_t) {
				sprintf(file_name1,"SLIR_annealR_%d.txt", temp_t);
				FILE *fp1 = fopen(file_name1, "a");
				fprintf(fp1, "%d\n", number_R);
				fclose(fp1);
				sprintf(file_name2,"SLIR_annealL_%d.txt", temp_t);
				FILE *fp2 = fopen(file_name2, "a");
				fprintf(fp2, "%d\n", number_L);
				fclose(fp2);
				sprintf(file_name3,"SLIR_annealI_%d.txt", temp_t);
				FILE *fp3 = fopen(file_name3, "a");
				fprintf(fp3, "%d\n", number_I);
				fclose(fp3);
				// printf("temp_t:%d  time: %f  R:%d \n", temp_t, t, number_R);
				temp_t += 1;
			}
			total_rate = number_I + number_L*k*number_I*a/n + (n - number_I - number_L - number_R)*k*number_I*a*(2-a)/n;
			t += 1/total_rate;
			if(state[ran_node] == 11) {
				state[ran_node] = 2;
				number_I -= 1;
				number_R += 1;
			}
			else if(state[ran_node] == 0) {
				ran_num = genrand_real3();
				if(ran_num < 2*(1 - a)/(2 - a)) {
					state[ran_node] = 1;
					number_L += 1;
				}
				else {
					state[ran_node] = 11;
					number_I += 1;
				}
			}
			else if(state[ran_node] == 1) {
				state[ran_node] = 11;
				number_I += 1;
				number_L -= 1;
			}
		}
		rate[ran_node] = 0;
		for(j = 0; j < k; j++) {
			label[neighbor[j]] = 0;
		}
	}
	// printf("number_R:%d time:%f\n", number_R, t);
	free(state);
	free(rate);
	free(neighbor);
	free(label);
	return number_R;
}

int SLIRStaticAveragetime(struct Link *net, int n, double a)
{
	int i, number_I, number_L, number_R, ran_node, temp, max_degree;
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate, total_rate, max_Itime, max_Ltime;
	int max_I, max_L, b;
	struct Link *p;
	double ran, t = 0.;
	char file_name[35];
	b = (int)(a*100);
	sprintf(file_name, "SLIR_static_RRG_avertime_%d.txt", b);
	FILE *fp = fopen(file_name, "a");
	max_degree = 6;
	max_I = 1;
	max_L = 0;
	max_Itime = 0.;
	max_Ltime = 0.;
	ran_node = (int)(n*genrand_real3());
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	number_L = 0;
	total_rate = 0;
	for(i = 0; i < n; i++) {
		if(state[i] == 11) rate[i] = 1;
		else {
			p = &net[i];
			temp = 0;
			while(p->next != NULL) {
				p = p->next;
				if(state[p->num] == 11) {
					temp += 1;
				}
			}
			rate[i] = temp*a*(2 - a);
		}
		total_rate += rate[i];
	}
	if(a*max_degree*(2 - a) > 1.0) max_rate = max_degree*a*(2 - a);
	else max_rate = 1.0;
	while(number_I > 0) {
		ran_node = (int)(n*genrand_real3());
		while(rate[ran_node] < 0.00001) {
			ran_node = (int)(n*genrand_real3());
		}
		// t += 1.0/total_rate; //73
		if(genrand_real3() < rate[ran_node]/max_rate) {
			t += 1.0/total_rate; //12
			if(number_I > max_I) {
				max_I = number_I;
				max_Itime = t;
			}
			if(number_L > max_L) {
				max_L = number_L;
				max_Ltime = t;
			}
			if(state[ran_node] == 11) {
				state[ran_node] = 2;
				total_rate -= rate[ran_node];
				rate[ran_node] = 0;
				number_I -= 1;
				number_R += 1;
				p = &net[ran_node];
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 0) {
						rate[p->num] -= a*(2 - a);
						total_rate -= a*(2 - a);
					}
					else if(state[p->num] == 1) {
						rate[p->num] -= a;
						total_rate -= a;
					} 
				}
			}
			else if(state[ran_node] == 0) {
				ran = genrand_real3();
				if(ran < 2*(1 - a)/(2 - a)) {
					state[ran_node] = 1;
					number_L += 1;
					p = &net[ran_node];
					temp = 0;
					while(p->next != NULL) {
						p = p->next;
						if(state[p->num] == 11) temp += 1;
					}
					rate[ran_node] = temp*a;
					total_rate += temp*a;
				}
				else {
					state[ran_node] = 11;
					rate[ran_node] = 1;
					total_rate += 1;
					number_I += 1;
					p = &net[ran_node];
					while(p->next != NULL) {
						p = p->next;
						if(state[p->num] == 0) {
							rate[p->num] += a*(2 - a);
							total_rate += a*(2 - a);
						} 
						else if(state[p->num] == 1) {
							rate[p->num] += a;
							total_rate += a;
						} 
					}
				}
			}
			else if(state[ran_node] == 1) {
				state[ran_node] = 11;
				rate[ran_node] = 1;
				total_rate += 1;
				number_I += 1;
				number_L -= 1;
				p = &net[ran_node];
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 0) {
						rate[p->num] += a*(2 - a);
						total_rate += a*(2 - a);
					} 
					else if(state[p->num] == 1) {
						rate[p->num] += a;
						total_rate += a;
					} 
				}
			}
		}
	}
	if(number_R > 1000) {
		fprintf(fp, "%f %f\n", max_Itime, max_Ltime);
	}
	// printf("number_R:%d \n", number_R);
	// printf("max_I:%d t:%f; max_L:%d, t:%f\n", max_I, max_Itime, max_L, max_Ltime);
	// printf("finaltime:%f\n", t);
	free(state);
	free(rate);
	fclose(fp);
	return number_R;
}

int SLIRAnnealAveragetime(int n, double a, int k)
{
	int i, j, number_I, number_L, number_R, ran_node;
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate, max_Itime, max_Ltime;
	int ran, k1, b, max_I, max_L;
	double ran_num;
	int *neighbor = (int *)malloc(k*sizeof(int));
	int *label = (int *)malloc(n*sizeof(int));
	double total_rate, t = 0.;
	char file_name[35];
	b = (int)(a*100);
	sprintf(file_name, "SLIR_anneal_RRG_avertime_%d.txt", b);
	FILE *fp = fopen(file_name, "a");
	max_I = 1;
	max_L = 0;
	max_Itime = 0.;
	max_Ltime = 0.;
	ran_node = (int)(n*genrand_real3());
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		rate[i] = 0;
		label[i] = 0;
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	number_L = 0;
	if(a*k*(2 - a) > 1.0) max_rate = k*a*(2 - a);
	else max_rate = 1.0;
	while(number_I > 0) {
		ran_node = (int)(n*genrand_real3());
		while(state[ran_node] == 2) {
			ran_node = (int)(n*genrand_real3());
		}
		if(state[ran_node] == 11) rate[ran_node] = 1;
		else if((state[ran_node] == 0) || (state[ran_node] == 1)) {
			k1 = 0;
			//choose neighbor
			for(j = 0; j < k; j++) {
				ran = (int)(n*genrand_real3());
				while((label[ran] == 1) || (ran == ran_node)) {
					ran = (int)(n*genrand_real3());
				}
				neighbor[j] = ran;
				label[ran] = 1;
				if(state[ran] == 11) k1 += 1;
			}
			if(state[ran_node] == 0) rate[ran_node] = k1*a*(2-a);
			else rate[ran_node] = k1*a;
		}
		if(genrand_real3() < rate[ran_node]/max_rate) {
			total_rate = number_I + (1.0/n)*number_L*k*number_I*a + (1.0/n)*(n - number_I - number_L - number_R)*k*number_I*a*(2-a);
			t += 1/total_rate;
			if(number_I > max_I) {
				max_I = number_I;
				max_Itime = t;
			}
			if(number_L > max_L) {
				max_L = number_L;
				max_Ltime = t;
			}
			if(state[ran_node] == 11) {
				state[ran_node] = 2;
				number_I -= 1;
				number_R += 1;
			}
			else if(state[ran_node] == 0) {
				ran_num = genrand_real3();
				if(ran_num < 2*(1 - a)/(2 - a)) {
					state[ran_node] = 1;
					number_L += 1;
				}
				else {
					state[ran_node] = 11;
					number_I += 1;
				}
			}
			else if(state[ran_node] == 1) {
				state[ran_node] = 11;
				number_I += 1;
				number_L -= 1;
			}
		}
		rate[ran_node] = 0;
		for(j = 0; j < k; j++) {
			label[neighbor[j]] = 0;
		}
	}
	if(number_R > 1000) {
		fprintf(fp, "%f %f\n", max_Itime, max_Ltime);
	}
	// printf("number_R:%d time:%f\n", number_R, t);
	// printf("max_I:%d t:%f; max_L:%d, t:%f\n", max_I, max_Itime, max_L, max_Ltime);
	free(state);
	free(rate);
	free(neighbor);
	free(label);
	fclose(fp);
	return number_R;
}

int main(int argc, char *argv[])
{
	int i, j, n, average_k = 6, number_R;
	double a;
	// struct Link *net;
	char file_name[25];
	int b;
	struct Link *p;
	n = 100000;
	// net = (struct Link *)malloc(n*sizeof(struct Link));
	if(argc > 1) {
		a = atof(argv[1]);
	} 
	else {
		printf("no parameter input\n");
		exit(1);
	}
	// b = (int)(a*100);
	// sprintf(file_name, "SLIR_anneal_RRG_%d.txt", b);
	// FILE *fp = fopen(file_name, "w");
	unsigned long idum = (unsigned long)time(NULL);
	init_genrand(idum);
	// for(i = 0; i < n; i++) {
	// 	net[i].num = i;
	// 	net[i].next = NULL;
	// }
	
	// for(i = 0; i < 100; i++) {
	// 	RegularRandomG(net, n, average_k);
	// 	for(j = 0; j < 100; j++) {
	// 		number_R = SLIRStaticLattice(net, n, a);
	// 		fprintf(fp, "%d\n", number_RL);
	// 	}
	// 	FreeNetwork(net, n);
	// }

	// for(i = 0; i < 10000; i++) {	
	// 	number_R = SLIRAnnealLattice(n, a, average_k);
	// 	fprintf(fp, "%d\n", number_R);
	// }
	// RegularRandomG(net, n, 6);
	for(i = 0; i < 10; i++) {
		// SLIRStaticAveragetime(net, n, a);
		SLIRAnnealAveragetime(n, a, 6);
	}
	// SLIRAnnealAveragetime(n, a, 6);
	// fclose(fp);
	return 0;
}