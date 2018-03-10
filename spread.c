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

void SaveNetwork(struct Link *net, int n)
{
    int i, temp;
    int *degree = (int *)malloc(n*sizeof(int));
    struct Link *p;
    FILE *fp = fopen("network_information1.txt", "w");
    for(i = 0; i < n; i++) {
        temp = 0;
        p = &net[i];
        while(p->next != NULL) {
            p = p->next;
            temp += 1;
        }
        degree[i] = temp;
    }
    for(i = 0; i < n; i++) {
        p = &net[i];
        fprintf(fp, "%d ", degree[i]);
        while(p->next != NULL) {
            p = p->next;
            fprintf(fp, "%d ", p->num);
        }
        fprintf(fp, "\n");
    }
    free(degree);
    fclose(fp);
}

void ReadNetwork(struct Link *net, int n)
{
    int i, temp, node;
    FILE *fp;
    if((fp = fopen("network_information1.txt","r")) == NULL){
        printf("cannot open this file.\n");
        exit(0);
    }
    for(i = 0; i < n; i++) {
        fscanf(fp, "%d", &temp);
        while(temp > 0) {
            fscanf(fp, "%d", &node);
            AddNode(&net[i], node);
            temp -= 1;
        }
    }
    fclose(fp);
}

//k is average degree of the network
void ERNetwork(struct Link *net, int n, int k)
{
	int i, node1, node2, total_link;
	struct Link *p;
	total_link = n*k/2;
	i = 0;
	while(i < total_link){
		node1 = (int)(n*genrand_real3());
		node2 = (int)(n*genrand_real3());
		while(node2 == node1) {
			node2 = (int)(n*genrand_real3());
		}
		p = &net[node1];
		while(p->next != NULL) {
			if(p->next->num == node2)  {
				break;
			}
			p = p->next;
		}
		if(p->next == NULL) {
			AddLink(net, node1, node2);
			i++;
		}
	}
}

//periodic boundary condition
void TriangularLattice(struct Link *net, int n) 
{
	int i, j, l;
	double size = (double)n;
	l = sqrt(size);
	for(i = 0; i < l; i++) {
		for(j = 0; j < l; j++) {
			AddLink(net, i*l + j, i*l + (j + 1)%l);
			AddLink(net, i*l + j, ((i + 1)%l)*l + j);
			AddLink(net, i*l + j, ((i + 1)%l)*l + (j + 1)%l);
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

long int fact(int n)
{
    int i;
    long int ans;
    ans = 1;
    for(i = n; i > 1; i--){
        ans = ans*i;
    }
    return ans;
}

//state: 'S'--'0', 'I'--'11', 'R'--'2'
//max rate
int SIRStatic(struct Link *net, int n, double a)
{
	int i, number_I, number_R, ran_node, temp, max_degree = 0;
	int *degree = (int *)malloc(n*sizeof(int));
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate = 0;
	struct Link *p;
	//calculate nodes' degree and choose a seed whose degree isn't zero
	for(i = 0; i < n; i++) {
		temp = 0;
		p = &net[i];
		while(p->next != NULL) {
			p = p->next;
			temp += 1;
		}
		degree[i] = temp;
		if(degree[i] > max_degree) max_degree = degree[i];
	}
	ran_node = (int)(n*genrand_real3());
	while(degree[ran_node] == 0) {
		ran_node = (int)(n*genrand_real3());
	}
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
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
			rate[i] = temp*a;
		}
	}
	if(a*max_degree > 1.0) max_rate = a*max_degree;
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
					if(state[p->num] == 0) rate[p->num] -= a;
				}
			}
			else if(state[ran_node] == 0) {
				state[ran_node] = 11;
				rate[ran_node] = 1;
				number_I += 1;
				p = &net[ran_node];
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 0) {
						rate[p->num] += a;
					}
				}
			}
			else printf("choose wrong\n");
		}
	}
	free(degree);
	free(state);
	free(rate);
	return number_R;
}

//state: 'S'--'0', 'I'--'11', 'R'--'2'
//max rate
int SIRAnneal(struct Link *net, int n, double a)
{
	int i, number_I, number_R, ran_node, temp, max_degree = 0, sum_deg, temp_neighbor;
	int *degree = (int *)malloc(n*sizeof(int));
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	int *index_neighbor = (int *)malloc(n*sizeof(int));
	int *neighbor;
	double ran, max_rate = 0;
	struct Link *p;
	//calculate nodes' degree and choose a seed whose degree isn't zero
	for(i = 0; i < n; i++) {
		temp = 0;
		p = &net[i];
		while(p->next != NULL) {
			p = p->next;
			temp += 1;
		}
		degree[i] = temp;
		if(degree[i] > max_degree) max_degree = degree[i];
	}
	neighbor = (int *)malloc(max_degree*sizeof(int));
	ran_node = (int)(n*genrand_real3());
	while(degree[ran_node] == 0) {
		ran_node = (int)(n*genrand_real3());
	}
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		index_neighbor[i] = 0;
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	for(i = 0; i < n; i++) {
		if(state[i] == 11) rate[i] = 1 + degree[i]*a;
		else rate[i] = 0.;
	}
	max_rate = 1 + a*max_degree;
	while(number_I > 0) {
		ran_node = (int)(n*genrand_real3());
		while(state[ran_node] != 11) {
			ran_node = (int)(n*genrand_real3());
		}
		if(genrand_real3() < rate[ran_node]/max_rate) {
			if(genrand_real3() < 1.0/rate[ran_node]) {
				state[ran_node] = 2;
				rate[ran_node] = 0;
				number_I -= 1;
				number_R += 1;
			}
			else {
				sum_deg = 0;
				for(i = 0; i < degree[ran_node]; i++) {
					temp_neighbor = (int)(n*genrand_real3());
					while((index_neighbor[temp_neighbor] == 1) || (temp_neighbor == ran_node)) {
						temp_neighbor = (int)(n*genrand_real3());
					}
					neighbor[i] = temp_neighbor;
					index_neighbor[temp_neighbor] = 1;
					sum_deg += degree[temp_neighbor] - 1;
				}
				ran = genrand_real3();
				for(i = 0; i < degree[ran_node]; i++) {
					if((degree[neighbor[i]] - 1)*1.0/sum_deg < ran) {
						ran -= (degree[neighbor[i]] - 1)*1.0/sum_deg;
					}
					else break;
				}
				if(state[neighbor[i]] == 0) {
					state[neighbor[i]] = 11;
					rate[neighbor[i]] = 1 + degree[neighbor[i]]*a;
					number_I += 1;
				}
				for(i = 0; i < degree[ran_node]; i++) {
					index_neighbor[neighbor[i]] = 0;
				}
			}
		}
	}
	free(degree);
	free(state);
	free(rate);
	free(index_neighbor);
	// printf("%d\n", number_R);
	return number_R;
}

//state: 'S'--'0', 'L'--'1', 'I'--'11', 'R'--'2'
int SLIRStaticER(struct Link *net, int n, double a)
{
	int i, number_I, number_L, number_R, ran_node, temp, max_degree = 0;
	int *degree = (int *)malloc(n*sizeof(int));
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate;
	struct Link *p;
	double ran;
	//calculate nodes' degree and choose a seed whose degree isn't zero
	for(i = 0; i < n; i++) {
		temp = 0;
		p = &net[i];
		while(p->next != NULL) {
			p = p->next;
			temp += 1;
		}
		degree[i] = temp;
		if(degree[i] > max_degree) max_degree = degree[i];
	}
	ran_node = (int)(n*genrand_real3());
	while(degree[ran_node] == 0) {
		ran_node = (int)(n*genrand_real3());
	}
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	number_L = 0;
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
					number_L += 1;
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
				number_L -= 1;
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
	free(degree);
	free(state);
	free(rate);
	return number_R;
}

//state: 'S'--'0', 'L'--'1', 'I'--'11', 'R'--'2'
//for regular graph such as lattice and RRG 
int SLIRStaticLattice(struct Link *net, int n, double a)
{
	int i, number_I, number_L, number_R, ran_node, temp, max_degree;
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
	number_L = 0;
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
					number_L += 1;
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
				number_L -= 1;
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
	return number_R;
}

//state: 'S'--'0','L'--'1', 'I'--'11', 'R'--'2'
//anneal RRG
int SLIRAnnealLattice(int n, double a, int k)
{
	int i, j, number_I, number_RL, number_R, ran_node;
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
	number_RL = 0;
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
				}
				else {
					state[ran_node] = 11;
					number_I += 1;
				}
			}
			else if(state[ran_node] == 1) {
				state[ran_node] = 11;
				number_I += 1;
				number_RL += 1;
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

int SLIRAnnealER(int n, double a, int k)
{
	int i, j, number_I, number_L, number_R, ran_node, temp, sum_deg, temp_neighbor;
	int *degree = (int *)malloc(n*sizeof(int));
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	int *index_neighbor = (int *)malloc(n*sizeof(int));
	int *neighbor;
	int *num = (int *)malloc((3*k + 1)*sizeof(int));
	double ran, max_rate, sum = 0;
	temp = 0;
	//calculate nodes' degree and choose a seed whose degree isn't zero
	for(i = 1; i <= 3*k; i++) {
		sum += pow(k, i)*exp(-k)/fact(i);
	}
	for(i = 1; i <= 3*k; i++) {
		num[i] = (int)(n*pow(k, i)*exp(-k)/(fact(i)*sum) + 0.5);
        temp += num[i];
	}
	num[1] += n - temp;
	temp = 0;
	for(i = 1; i <= 3*k; i++) {
		for(j = 0; j < num[i]; j++) {
			degree[temp] = i;
			temp += 1;
		}
	}
	// neighbor = (int *)malloc(max_degree*sizeof(int));
	ran_node = (int)(n*genrand_real3());
	while(degree[ran_node] == 0) {
		ran_node = (int)(n*genrand_real3());
	}
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		index_neighbor[i] = 0;
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_L = 0;
	number_R = 0;
	for(i = 0; i < n; i++) {
		if(state[i] == 11) rate[i] = 1 + degree[i]*a*(2-a);
		else rate[i] = 0.;
	}
	max_rate = 1 + 3*k*a*(2-a);
	while(number_I > 0) {
		ran_node = (int)(n*genrand_real3());
		while(state[ran_node] != 11) {
			ran_node = (int)(n*genrand_real3());
		}
		if(genrand_real3() < rate[ran_node]/max_rate) {
			if(genrand_real3() < 1.0/rate[ran_node]) {
				state[ran_node] = 2;
				rate[ran_node] = 0;
				number_I -= 1;
				number_R += 1;
			}
			else {
				sum_deg = 0;
				neighbor = (int *)malloc(degree[ran_node]*sizeof(int));
				for(i = 0; i < degree[ran_node]; i++) {
					temp_neighbor = (int)(n*genrand_real3());
					while((index_neighbor[temp_neighbor] == 1) || (temp_neighbor == ran_node)) {
						temp_neighbor = (int)(n*genrand_real3());
					}
					neighbor[i] = temp_neighbor;
					index_neighbor[temp_neighbor] = 1;
					sum_deg += degree[temp_neighbor] - 1;
				}
				ran = genrand_real3();
				for(i = 0; i < degree[ran_node]; i++) {
					if((degree[neighbor[i]] - 1)*1.0/sum_deg < ran) {
						ran -= (degree[neighbor[i]] - 1)*1.0/sum_deg;
					}
					else break;
				}
				if(state[neighbor[i]] == 1) {
					if(genrand_real3() < 1/(2-a)) {
						state[neighbor[i]] = 11;
						rate[neighbor[i]] = 1 + degree[neighbor[i]]*a*(2-a);
						number_I += 1;
						number_L -= 1;
					}
				}
				else if(state[neighbor[i]] == 0){
					if(genrand_real3() < a/(2-a)) {
						state[neighbor[i]] = 11;
						rate[neighbor[i]] = 1 + degree[neighbor[i]]*a*(2-a);
						number_I += 1;
					}
					else {
						state[neighbor[i]] = 1;
						number_L += 1;
					}
				}
				for(i = 0; i < degree[ran_node]; i++) {
					index_neighbor[neighbor[i]] = 0;
				}
				free(neighbor);
			}
		}
	}
	free(num);
	free(degree);
	free(state);
	free(rate);
	free(index_neighbor);
	// printf("%d\n", number_R);
	return number_R;
}

//state: 'S'--'0','L'--'1', 'I'--'11', 'R'--'2'
int SLIRAnnealERmy(int n, double a, int k) 
{
	int i, j, number_L, number_I, number_R, ran_node, temp = 0, total_deg = 0;
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate;
	int ran, k1;
	double sum = 0;
	int *num = (int *)malloc((3*k + 1)*sizeof(int));
	int *neighbor;
	int *label = (int *)malloc(n*sizeof(int));
	int *degree = (int *)malloc(n*sizeof(int));
	int *choose_neighbor;
	//assign degree to each node and choose a seed node
	for(i = 1; i <= 3*k; i++) {
		sum += pow(k, i)*exp(-k)/fact(i);
	}
	for(i = 1; i <= 3*k; i++) {
		num[i] = (int)(n*pow(k, i)*exp(-k)/(fact(i)*sum) + 0.5);
        temp += num[i];
	}
	num[1] += n - temp;
	temp = 0;
	for(i = 1; i <= 3*k; i++) {
		for(j = 0; j < num[i]; j++) {
			degree[temp] = i;
			temp += 1;
		}
		total_deg += i*num[i];
	}
	choose_neighbor = (int *)malloc(total_deg*sizeof(int));
	temp = 0;
	for(i = 0; i < n; i++) {
		for(j = 0; j < degree[i]; j++) {
			choose_neighbor[temp] = i;
			temp += 1;
		}
	}
	ran = (int)(n*genrand_real3());
	//initialize nodes' state and start to spread
	for(i = 0; i < n; i++) {
		rate[i] = 0;
		label[i] = 0;
		if(i == ran) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_L = 0;
	number_R = 0;
	if(3*k*a*(2 - a) > 1.0) max_rate = 3*k*a*(2 - a);
	else max_rate = 1.0;
	while(number_I > 0) {
		ran_node = (int)(n*genrand_real3());
		while(state[ran_node] == 2) {
			ran_node = (int)(n*genrand_real3());
		}
		if(state[ran_node] == 11) rate[ran_node] = 1;
		else if((state[ran_node] == 0) || (state[ran_node] == 1)) {
			k1 = 0;
			neighbor = (int *)malloc(degree[ran_node]*sizeof(int));
			// choose neighbor
			for(j = 0; j < degree[ran_node]; j++) {
				ran = (int)(total_deg*genrand_real3());
				ran = choose_neighbor[ran];
				while((label[ran] == 1) || (ran == ran_node)) {
					ran = (int)(total_deg*genrand_real3());
					ran = choose_neighbor[ran];
				}
				neighbor[j] = ran;
				label[ran] = 1;
				if(state[ran] == 11) k1 += 1;
			}
			if(state[ran_node] == 0) rate[ran_node] = k1*a*(2-a);
			else rate[ran_node] = k1*a;
			for(j = 0; j < degree[ran_node]; j++) {
				label[neighbor[j]] = 0;
			}
			free(neighbor);
		}
		if(genrand_real3() < rate[ran_node]/max_rate) {
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
	}
	// printf("%d\n", number_R);
	free(num);
	free(state);
	free(rate);
	free(degree);
	free(label);
	printf("%d\n", number_R);
	// fclose(fp);
	return number_R;
}

//state: 'S'--'0', 'L'--'1', 'I'--'11', 'R'--'2'
int SLIRwenxian1(struct Link *net, int n, double a)
{
	int i, number_I, number_L, number_R, ran_node, temp, max_degree = 0;
	int *degree = (int *)malloc(n*sizeof(int));
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate = 0;
	struct Link *p;
	//calculate nodes' degree and choose a seed whose degree isn't zero
	for(i = 0; i < n; i++) {
		temp = 0;
		p = &net[i];
		while(p->next != NULL) {
			p = p->next;
			temp += 1;
		}
		degree[i] = temp;
		if(degree[i] > max_degree) max_degree = degree[i];
	}
	ran_node = (int)(n*genrand_real3());
	while(degree[ran_node] == 0) {
		ran_node = (int)(n*genrand_real3());
	}
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	number_L = 0;
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
			rate[i] = temp*a;
		}
	}
	if(a*max_degree > 1.0) max_rate = max_degree*a;
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
					if(state[p->num] == 0) rate[p->num] -= a;
					else if(state[p->num == 1]) rate[p->num] -= 1;
				}
			}
			else if(state[ran_node] == 0) {
				state[ran_node] = 1;
				number_L += 1;
				p = &net[ran_node];
				temp = 0;
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 11) temp += 1;
				}
				rate[ran_node] = temp;
			}
			else if(state[ran_node] == 1) {
				state[ran_node] = 11;
				number_I += 1;
				number_L -= 1;
				p = &net[ran_node];
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 0) rate[p->num] += a;
					else if(state[p->num == 1]) rate[p->num] += 1;
				}
			}
		}
	}
	// printf("number_R:%d\n", number_R);
	free(degree);
	free(state);
	free(rate);
	return number_R;
}

//state: 'S'--'0', 'L'--'1', 'I'--'11', 'R'--'2'
int SLIRwenxian2(struct Link *net, int n, double a)
{
	int i, number_I, number_L, number_R, ran_node, temp, max_degree = 0;
	int *degree = (int *)malloc(n*sizeof(int));
	int *state = (int *)malloc(n*sizeof(int));
	double *rate = (double *)malloc(n*sizeof(double));
	double max_rate = 0;
	struct Link *p;
	//calculate nodes' degree and choose a seed whose degree isn't zero
	for(i = 0; i < n; i++) {
		temp = 0;
		p = &net[i];
		while(p->next != NULL) {
			p = p->next;
			temp += 1;
		}
		degree[i] = temp;
		if(degree[i] > max_degree) max_degree = degree[i];
	}
	ran_node = (int)(n*genrand_real3());
	while(degree[ran_node] == 0) {
		ran_node = (int)(n*genrand_real3());
	}
	//initialize nodes' state and build the rate list
	for(i = 0; i < n; i++) {
		if(i == ran_node) state[i] = 11;
		else state[i] = 0;
	}
	number_I = 1;
	number_R = 0;
	number_L = 0;
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
			rate[i] = temp*a;
		}
	}
	if(a*max_degree > 1.0) max_rate = max_degree*a;
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
					if((state[p->num] == 0) || (state[p->num] == 1)) rate[p->num] -= a;
				}
			}
			else if(state[ran_node] == 0) {
				state[ran_node] = 1;
				number_L += 1;
				p = &net[ran_node];
				temp = 0;
				while(p->next != NULL) {
					p = p->next;
					if(state[p->num] == 11) temp += 1;
				}
				rate[ran_node] = temp*a;
			}
			else if(state[ran_node] == 1) {
				state[ran_node] = 11;
				rate[ran_node] = 1;
				number_I += 1;
				number_L -= 1;
				p = &net[ran_node];
				while(p->next != NULL) {
					p = p->next;
					if((state[p->num] == 0) || (state[p->num] == 1)) rate[p->num] += a;
				}
			}
		}
	}
	// printf("number_R:%d\n", number_R);
	free(state);
	free(rate);
	free(degree);
	return number_R;
}

int main(int argc, char *argv[])
{
	int i, j, n = 100000, average_k = 6, number_R;
	double a;
	struct Link *net = (struct Link *)malloc(n*sizeof(struct Link));
	char file_name[25];
	int b;
	struct Link *p;
	if(argc > 1) {
		a = atof(argv[1]);
	} 
	else {
		printf("no parameter input\n");
		exit(1);
	}
	// b = (int)(a*100);
	// sprintf(file_name, "SLIR_anneal_ER_%d.txt", b);
	// FILE *fp = fopen(file_name, "w");
	unsigned long idum = (unsigned long)time(NULL);
	init_genrand(idum);
	// for(i = 0; i < n; i++) {
	// 	net[i].num = i;
	// 	net[i].next = NULL;
	// }
	// for(i = 0; i < 100; i++) {
	// 	ERNetwork(net, n, average_k);
	// 	for(j = 0; j < 100; j++) {
	// 		number_R = SLIRannealER(net, n, a);
	// 		fprintf(fp, "%d\n", number_R);
	// 	}
	// 	FreeNetwork(net, n);
	// }
	// TriangularLattice(net, n);
	// for(i = 0; i < 10000; i++) {	
	// 	number_R = SLIRStaticLattice(net, n, a);
	// 	fprintf(fp, "%d\n", number_R);
	// }
	// FreeNetwork(net, n);
	// ERNetwork(net, n, average_k);
	// SLIRAnnealERtry(n, a, 6);
	for(i = 0; i < 100; i++) {
		number_R = SLIRAnnealERtry(n, a, 6);
		if(number_R > 1000) break;
	}
	
	// fclose(fp);
	return 0;
}