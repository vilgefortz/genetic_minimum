//#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <stropts.h>

using namespace std;

struct context
{
    float* values;
    float* stack;
    float* tmp;
    int pointer;
    int tmp_pointer;
};

//operacja z listy wygenerowana z podanego wzoru
struct operand
{
    char* name;
    void (*cmd)(context&);
};



context ctx;
bool finished = false;
operand* op;
int plen = 0;
//potrzebne do bitowego przesuwania floata
typedef union
{
    long i;
    float f;
} u;

struct pv
{
    float* val; //genotyp
    float fit; //przystoowanie
};

struct para //repreze ntuje pare rodzicow
{
    float* a;
    float* b;
};

float target_time = 1000; //ms
int ARGS = 1;

int POPLEN = 0;
float* LOW;
float* HIGH;
void find_min();

float** parents;
float** children;
float** tmpparents;

float calculate(float*);
float** create_initial();
pv** calculate_fitness(float**, int);
float** get_parents(float**, pv**);
float** get_children(float**);
float* normalized = (float*)malloc (sizeof *normalized * POPLEN);
char* main_buffer;
char* buffer;
char* line;
char** parts;
size_t size = 999;

int t_start=0;
int t_len=0;

char* sub_string (const char* input, int offset, int len, char* dest)
{
    int input_len = strlen (input);

    if (offset + len > input_len)
    {
        return NULL;
    }

    strncpy (dest, input + offset, len);
    dest[len] = '\0';
    return dest;
}

bool read_word(char* line)
{
    while (line[t_start] != '\n' && line[t_start] != '\r')
    {
        if (line[t_start] == ' ' || line[t_start] == '\t')
        {
            t_start++;
            continue;
        }
        else
        {
            break;
        }
    }
    if (line[t_start] == '\n' || line[t_start] == '\0' || line[t_start] == '\r')
    {
        t_start = 0;
        t_len = 0;
        return false;
    }
    while (line[t_start + t_len] != ' ' && line[t_start + t_len] != '\n' && line[t_start + t_len] != '\t' && line[t_start] != '\0' && line[t_start + t_len] != '\r')
    {
        t_len++;
    }
    sub_string(line, t_start, t_len, buffer);
    t_start = t_start + t_len;
    t_len = 0;
    return true;
}

void get_line() {
    t_start=0;
    t_len=0;
    do {
        getline(&line, &size, stdin);
    } while (line[0] == '#');
}

void read_input()
{
    buffer = (char*)malloc(sizeof *buffer * 1000);
    line = (char*)malloc(sizeof *line * 50000);
    parts = (char**)malloc(sizeof **parts * 300);
    for (int i = 0; i< 300; i++)
    {
        parts[i] = (char*)malloc(sizeof *parts * 1000);
    }
    get_line();
    ////cout << line;
    read_word(line);
    ARGS = atoi(buffer);
    ////cout << ARGS << endl;
    read_word(line);
    target_time = atoi(buffer);
    LOW = (float*) malloc(sizeof *LOW * ARGS);
    HIGH = (float*) malloc(sizeof *LOW * ARGS);
    for (int i = 0; i < ARGS; i++)
    {
        get_line();
        ////cout << line;
        read_word(line);
        LOW[i] = atof(buffer);
        ////cout << "LOW " << buffer << " = "<< LOW[i] <<endl;
        read_word(line);
        HIGH[i] = atof(buffer);
        ////cout << "HIGH "<< buffer << " = " << HIGH[i]<<endl;
    }
    plen = 0;
    get_line();
    ////cout << line;
    bool more = true;
    while (more && read_word(line)) {
        more = false;
    do
    {
        ////cout << " PARTS : " << buffer << " " << strlen(buffer) << endl;
        for (int i=0 ; i < strlen(buffer); i++) {
            ////cout << "SIGN : " << (int)buffer[i];
        }
        strcpy(parts[plen++], buffer);
    } while (read_word(line));
    if (!feof(stdin)) {
        get_line();
        more = true;
    }
    }
    //exit (1);
}

float space(float* x)
{
    return calculate(x);
}

float norm(float i, float low, float high)
{
    return low + (high-low) * i;
}

float* norm_all(float* v)
{
    for (int i =0 ; i< ARGS; i++)
    {
        normalized[i] = norm(v[i], LOW[i], HIGH[i]);
    }
    return normalized;
}

float space_norm(float* i)
{
    return space(norm_all(i));
}
void polish();
int main()
{
    read_input();
    polish();
    //exit(0);
    for (int i=0; i< ARGS; i++)
    {
        POPLEN += HIGH[i] - LOW[i];
    }
    POPLEN *= ARGS;
    //if (POPLEN < 100) POPLEN = 100;
    //if (POPLEN > 2000)
    POPLEN = 2000;
    tmpparents = (float**) malloc (sizeof *tmpparents * POPLEN) ;
    parents = (float**) malloc (sizeof *parents * POPLEN) ;
    children = (float**) malloc (sizeof *children * POPLEN) ;
    srand(time(NULL));
    for (int i=0; i<POPLEN; i++)
    {
        tmpparents[i] = (float*) malloc (sizeof **tmpparents * ARGS) ;
        parents[i] = (float*) malloc (sizeof **parents * ARGS) ;
        children[i] = (float*) malloc (sizeof **children * ARGS) ;
    }
    find_min();

}

void find_min()
{
    float min = FLT_MAX;
    float val = 0;
    float** pop = create_initial();
    float elapsed = 0;
    clock_t t1,t2;
    float diff;
    int iter = 0;
    do
    {

        t1 = clock();
        pv** fitness = calculate_fitness(pop, POPLEN);
        ////////cout << "getting parents" << endl;
        float** parents = get_parents(pop, fitness);
        ////////cout << "getting children" << endl;
        pop = get_children(parents);
        t2 = clock();
        diff = ((float)(t2 - t1) / (CLOCKS_PER_SEC) );
        elapsed += diff;
        iter ++;
        //////cout << "ITER " << iter << endl;
    }
    while (!finished && (target_time - diff) > elapsed);
    float* res;
    for (int i = 0; i<POPLEN; i++)
    {
        if (space_norm(pop[i]) < min)
        {
            min = space_norm(pop[i]);
            res = pop[i];
        }
    }
    for (int i=0; i<ARGS; i++)
    {
        printf ("%f ", norm(res[i], LOW[i], HIGH[i]));
    }
    printf("\n%d", iter);

}

int sgn (float x)
{
    return (x > 0) - (x < 0);
}

int cmp_pv (const void* l, const void* r)
{
    return sgn((*(pv**)r)->fit - (*(pv**)l)->fit);
}

float** get_parents(float** pop, pv** fitness)
{
    //TODO czy potrzebne?
    //float** parents = new float[POPLEN];

    //uzycia parentow na poszcegolnych indeksach
    //int uses[POPLEN];
    float fitnessD[POPLEN];
    pv** value_fitness;
    value_fitness = (pv**)malloc (sizeof *value_fitness * POPLEN);
    //uses[0] = 0;
    fitnessD[0] = fitness[0]->fit;
    value_fitness[0] = (pv*)malloc(sizeof **value_fitness);

    //przygotowanie zmiennych
    for (int i = 1; i < POPLEN; i++)
    {
        value_fitness[i] = (pv*)malloc(sizeof **value_fitness);
        //    value_fitness[i] = new pv;
       // uses[i] = 0;
        fitnessD[i] = fitness[i]->fit + fitnessD[i - 1];
    }
    float sum = fitnessD[POPLEN-1];
    for (int i = 0; i < POPLEN; i++)
    {
            float roulette = (float)rand()/RAND_MAX * sum;
           // found = false;
            value_fitness[i]->val = pop[0];
            value_fitness[i]->fit = fitness[0]->fit;
            for (int j = 1; j < POPLEN; j++)
            {
                if (fitnessD[j - 1] <= roulette && roulette < fitnessD[j])
                {
                    value_fitness[i]->val = pop[j];
                    value_fitness[i]->fit = fitness[j]->fit;
                    break;
                }
            }
    }

    qsort(value_fitness, POPLEN ,sizeof value_fitness, cmp_pv);
    for (int i = 0; i < POPLEN; i++)
    {
        parents[i] = value_fitness[i]->val;
    }
    return parents;
}


float** create_initial()
{
    float** pop = parents;
    for (int k = 0 ; k < ARGS; k++)
    {
        pop[0][k] = 0;
        pop[1][k] = 1;
        for (int i=2; i<POPLEN; i++)
        {
            pop[i][k] = (float)rand() / RAND_MAX;
            ////cout << pop[i][k] << "\n";
        }
    }

    return pop;
}

pv** calculate_fitness(float** pop, int len)
{
    pv** fitness;
    fitness = (pv**)malloc (sizeof *fitness * len);
    float abssum = 0;
    float sum = 0;
    for (int i = 0; i < len; i++)
    {
        fitness[i] = (pv*)malloc(sizeof **fitness);
        abssum += fabsf(space_norm(pop[i]));
        sum += space_norm(pop[i]);
    }
    for (int i = 0; i < len; i++)
    {
        float a = abssum - (float) space_norm (pop[i]);
        float b = sum + abssum * len;
        fitness[i]->val = pop[i];
        fitness[i]->fit = (a / b);
    }
    return fitness;
}

int mutate (int f)
{
    return (float)rand()/RAND_MAX < 0.2? f xor (1 << rand() % 23) : f;
}

float** get_children(float** parents)
{
    para* pairs = (para*)malloc(sizeof *pairs * POPLEN/2);
    if (pairs == NULL) {
        //cout << "MALLOC ERROR";
        exit(0);
    }
    int c = 0;
    for (int i = 0; i < POPLEN - 1; i++)
    {
        int uses = 0;
        if (c>=(POPLEN/2)) break;
        for (int j = i + 1; j < POPLEN; j++)
        {
            if (c>=(POPLEN/2) || parents[i] == parents[j] || (float)rand()/RAND_MAX < 0.5)
                continue;
            uses++;
            if (uses > 2)
            {
                uses = 0;
                i++;
                continue;
            }
            pairs[c].a = parents[i];
            pairs[c].b = parents[j];
            c++;
        }
    }
    if (c < POPLEN/2) {
        finished = true;
        return parents;
    }
    float** cvals = children;
    //pv* ch = new pv[c];
    ////////cout << "SIIIZE " << c << "\n";
    for (int k = 0; k< ARGS; k++)
    {
        for (int i = 0; i < POPLEN/2; i++)
        {
            //  //////cout << i << " "<< pairs[i].a;
            // //////cout << " " << pairs[i].b << "\n";
back:
            u maska,maskb,x,y,q,p;
            int width = rand() % 32;
            int mask = -1 << width;
            int imask = ~mask;
            maska.f = rand() / RAND_MAX;
            maskb.f = rand() / RAND_MAX;
            x.f = pairs[i].a[k];
            y.f = pairs[i].b[k];
            p.i = x.i & mask | y.i & imask;
            q.i = x.i & imask | y.i & mask;

            p.i = mutate(p.i);
            q.i = mutate(q.i);
            if (p.f > 1 || p.f < 0 || q.f > 1 || q.f < 0) goto back;
            cvals[i*2][k] = p.f;
            cvals[i*2 + 1][k] = q.f;
        }
    }
    free(pairs);
    return cvals;
}


//x, y, z, a, b, c, d, e, f, g

void push (context &ctx)
{
    ctx.stack[ctx.pointer++] = ctx.tmp[ctx.tmp_pointer++];
}

void push_x (context &ctx)
{
    ctx.stack[ctx.pointer++] = ctx.values[0];
}

void push_y (context &ctx)
{
    ctx.stack[ctx.pointer++] = ctx.values[1];
}

void push_z (context &ctx)
{
    ctx.stack[ctx.pointer++] = ctx.values[2];
}

void push_a (context &ctx)
{
    ctx.stack[ctx.pointer++] = ctx.values[3];
}

void push_b (context &ctx)
{
    ctx.stack[ctx.pointer++] = ctx.values[4];
}

void push_c (context &ctx)
{
    ctx.stack[ctx.pointer++] = ctx.values[5];
}

void push_d (context &ctx)
{
    ctx.stack[ctx.pointer++] = ctx.values[6];
}

void push_e (context &ctx)
{
    ctx.stack[ctx.pointer++] = ctx.values[7];
}

void push_f (context &ctx)
{
    ctx.stack[ctx.pointer++] = ctx.values[8];
}

void push_g (context &ctx)
{
    ctx.stack[ctx.pointer++] = ctx.values[9];
}

void add_a_b (context &ctx)
{
    ctx.stack[ctx.pointer-2] = ctx.stack[ctx.pointer-2] + ctx.stack[ctx.pointer---1];
}

void sub_a_b (context &ctx)
{
    ctx.stack[ctx.pointer-2] = ctx.stack[ctx.pointer-2] - ctx.stack[ctx.pointer---1];
}

void mul_a_b (context &ctx)
{
    ctx.stack[ctx.pointer-2] = ctx.stack[ctx.pointer-2] * ctx.stack[ctx.pointer---1];
}

void div_a_b (context &ctx)
{
    ctx.stack[ctx.pointer-2] = ctx.stack[ctx.pointer-2] / ctx.stack[ctx.pointer---1];
}

void pow_a (context &ctx)
{
    ctx.stack[ctx.pointer-2] = powf(ctx.stack[ctx.pointer-2], ctx.stack[ctx.pointer---1]);
}

void sin_a (context &ctx)
{
    ctx.stack[ctx.pointer-1] = sin(ctx.stack[ctx.pointer-1]);
}

void cos_a (context &ctx)
{
    ctx.stack[ctx.pointer-1] = cos(ctx.stack[ctx.pointer-1]);
}

void push_pi (context &ctx)
{
    ctx.stack[ctx.pointer++] = M_PI;
}

void neg_a (context &ctx)
{
    ctx.stack[ctx.pointer-1] = -(ctx.stack[ctx.pointer-1]);
}

void abs_a (context &ctx)
{
    ////cout << ctx.stack[ctx.pointer-1] << " ";
    ctx.stack[ctx.pointer-1] = fabsf(ctx.stack[ctx.pointer-1]);
    ////cout << fabsf(ctx.stack[ctx.pointer-1]) << " ";
}

void sqr_a (context &ctx)
{
    ////////cout << "POWERING VALUE " << ctx.stack[ctx.pointer-1] << endl;
    ctx.stack[ctx.pointer-1] = (ctx.stack[ctx.pointer-1]) * (ctx.stack[ctx.pointer-1]);
}

void sqrt_a (context &ctx)
{
    ctx.stack[ctx.pointer-1] = sqrtf(ctx.stack[ctx.pointer-1]);
}

void exp_a (context &ctx)
{
    ctx.stack[ctx.pointer-1] = expf(ctx.stack[ctx.pointer-1]);
}

void log_a (context &ctx)
{
    ctx.stack[ctx.pointer-1] = logf(ctx.stack[ctx.pointer-1]);
}




operand operands [25] =
{
    {
        "",
        push
    },

    {
        "x",
        push_x
    },
    {
        "y",
        push_y
    },
    {
        "z",
        push_z
    },
    {
        "a",
        push_a
    },
    {
        "b",
        push_b
    },
    {
        "c",
        push_c
    },
    {
        "d",
        push_d
    },
    {
        "e",
        push_e
    },
    {
        "f",
        push_f
    },
    {
        "g",
        push_g
    },
    {
        "+",
        add_a_b
    },
    {
        "-",
        sub_a_b
    },
    {
        "*",
        mul_a_b
    },
    {
        "/",
        div_a_b
    },
    {
        "pi",
        push_pi
    },
    {
        "neg",
        neg_a
    },
    {
        "abs",
        abs_a
    },
    {
        "sqr",
        sqr_a
    },
    {
        "sqrt",
        sqrt_a
    },
    {
        "exp",
        exp_a
    },
    {
        "log",
        log_a
    },
    {
        "pow",
        pow_a
    },
    {
        "sin",
        sin_a
    },
    {
        "cos",
        cos_a
    }

};


int olen = 25;
void polish()
{
//    context ctx;
    int tmp_pointer = 0;
    ctx.tmp_pointer = 0;
    ctx.pointer = 0;
    ctx.values = (float*) malloc(sizeof *(ctx.values) * 10);
    ctx.values[0] = 5;
    ctx.values[1] = 6;
    //plen = 9;
    ctx.stack = (float*) malloc(sizeof *(ctx.stack) * 200);
    ctx.tmp = (float*) malloc(sizeof *(ctx.tmp) * 200);
    op = (operand*)malloc(sizeof(*op) * plen);
    for (int i = 0; i< plen; i++)
    {
        bool found = false;
        for (int j = 0; j< olen; j++)    //operands amount
        {
            //    //////cout << operands[j].name << endl;

            if (strcmp(operands[j].name, parts[i]) == 0)
            {
                op[i] = operands[j];
                //cout << op[i].name << endl << flush;
                found = true;
                break;
            }
        }
        if (!found)
        {
            //////cout << "KONWERSJA" << parts[i] << " = " << atof(parts[i]) << endl;
            ctx.tmp[tmp_pointer++] = atof(parts[i]);
            //cout << "NUM :" << atof(parts[i]) << endl;
            op[i] = operands[0];
            //cout << op[i].name << endl  << flush;
        }
    }

    for (int i=0; i<plen; i++)
    {
        op[i].cmd(ctx);
        //////cout << op[i].name << " "<< ctx.stack[0] << " " <<  ctx.stack[1] << "\n";
    }
}

float calculate(float* v)
{
    ctx.pointer = 0;
    ctx.tmp_pointer = 0;
    ctx.values = v;
    for (int i=0; i<plen; i++)
    {
        op[i].cmd(ctx);
    }
    //  //////cout << ctx.stack[0] << endl;
    return ctx.stack[0];
}




