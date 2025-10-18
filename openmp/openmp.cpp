#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <vector>
#include <chrono>
#include <iostream>


int number_bacteria;
char** bacteria_name;
long M_6, M_5, M_4;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};
#define encode(ch)		code[ch-'A']
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010

void Init()
{
	M_4 = 1;
	for (int i=0; i<LEN-2; i++)	// M_4 = AA_NUMBER ^ (LEN-2);
		M_4 *= AA_NUMBER; 
	M_5 = M_4 * AA_NUMBER;		// M_5 = AA_NUMBER ^ (LEN-1);
	M_6  = M_5 *AA_NUMBER;			// M_6  = AA_NUMBER ^ (LEN);
}

class Bacteria
{
private:
    long *mers_6;
    long *mers_5;
    long mers_1[AA_NUMBER];
    long indexs;
    long total_mers_6;
    long total_mers_1;
    long complement;

    void InitVectors()
    {
        mers_6 = new long[M_6];
        mers_5 = new long[M_5];
        memset(mers_6, 0, M_6 * sizeof(long));
        memset(mers_5, 0, M_5 * sizeof(long));
        memset(mers_1, 0, AA_NUMBER * sizeof(long));
        total_mers_6 = 0;
        total_mers_1 = 0;
        complement = 0;
    }

    void init_buffer(char *buffer)
    {
        complement++;
        indexs = 0;
        for (int i = 0; i < LEN - 1; i++)
        {
            short enc = encode(buffer[i]);
            mers_1[enc]++;
            total_mers_1++;
            indexs = indexs * AA_NUMBER + enc;
        }
        mers_5[indexs]++;
    }

    void cont_buffer(char ch)
    {
        short enc = encode(ch);
        mers_1[enc]++;
        total_mers_1++;
        long index = indexs * AA_NUMBER + enc;
        mers_6[index]++;
        total_mers_6++;
        indexs = (indexs % M_4) * AA_NUMBER + enc;
        mers_5[indexs]++;
    }

public:
    long count;
    // double *sig_t_vec;
    // long *sig_t_indx_vec;

    std::vector<double> sig_t_vec;
    std::vector<long> sig_t_indx_vec;

    Bacteria(char *filename)
    {
        FILE *bacteria_file = fopen(filename, "r");
        if (bacteria_file == NULL)
        {
            fprintf(stderr, "Error: failed to open file %s\n", filename);
            exit(1);
        }

        InitVectors();

        char ch;
        while ((ch = fgetc(bacteria_file)) != EOF)
        {
            if (ch == '>')
            {
                while (fgetc(bacteria_file) != '\n')
                    ; // skip rest of line

                char buffer[LEN - 1];
                fread(buffer, sizeof(char), LEN - 1, bacteria_file);
                init_buffer(buffer);
            }
            else if (ch != '\n' && ch != '\r')
                cont_buffer(ch);
        }

        long total_plus_complement = total_mers_6 + complement;
        double total_div_2 = total_mers_6 * 0.5;
        int i_mod_aa_number = 0;
        int i_div_aa_number = 0;
        long i_mod_M1 = 0;
        long i_div_M1 = 0;

        double one_l_div_total[AA_NUMBER];
        for (int i = 0; i < AA_NUMBER; i++)
            one_l_div_total[i] = (double)mers_1[i] / total_mers_1;

        double *mers_5_div_total = new double[M_5];
        for (int i = 0; i < M_5; i++)
            mers_5_div_total[i] = (double)mers_5[i] / total_plus_complement;

        count = 0;

        // std::map<long, double> t;

        // sig_t_vec = new double[count];
        // sig_t_indx_vec = new long[count];

        for (long i = 0; i < M_6; i++)
        {
            double p1 = mers_5_div_total[i_div_aa_number];
            double p2 = one_l_div_total[i_mod_aa_number];
            double p3 = mers_5_div_total[i_mod_M1];
            double p4 = one_l_div_total[i_div_M1];
            double stochastic = (p1 * p2 + p3 * p4) * total_div_2;

            if (i_mod_aa_number == AA_NUMBER - 1)
            {
                i_mod_aa_number = 0;
                i_div_aa_number++;
            }
            else
                i_mod_aa_number++;

            if (i_mod_M1 == M_5 - 1)
            {
                i_mod_M1 = 0;
                i_div_M1++;
            }
            else
                i_mod_M1++;

            if (stochastic > EPSILON)
            {
                sig_t_vec.push_back((mers_6[i] - stochastic) / stochastic);
                sig_t_indx_vec.push_back(i);
                count++;
            }
        }

        delete mers_6;
        delete mers_5;

        fclose(bacteria_file);
    }
};

void ReadInputFile(const char* input_name)
{
	FILE* input_file = fopen(input_name, "r");
	if (input_file == NULL) {
		fprintf(stderr, "Error: failed to open file %s\n", input_name);
		exit(1);
	}

	fscanf(input_file, "%d", &number_bacteria);
	bacteria_name = new char*[number_bacteria];

	for(long i=0;i<number_bacteria;i++)
	{
		char name[10];
		fscanf(input_file, "%s", name);
		bacteria_name[i] = new char[30];
		snprintf(bacteria_name[i], 30, "../data/%s.faa", name);
	}
	fclose(input_file);
}

double CompareBacteria(Bacteria *b1, Bacteria *b2)
{
    double correlation = 0;
    double vector_len1 = 0;
    double vector_len2 = 0;
    long p1 = 0;
    long p2 = 0;
    while (p1 < b1->sig_t_vec.size() && p2 < b2->sig_t_vec.size())
    {
        long n1 = b1->sig_t_indx_vec[p1];
        long n2 = b2->sig_t_indx_vec[p2];
        if (n1 < n2)
        {
            double t1 = b1->sig_t_vec[p1];
            vector_len1 += (t1 * t1);
            p1++;
        }
        else if (n2 < n1)
        {
            double t2 = b2->sig_t_vec[p2];
            p2++;
            vector_len2 += (t2 * t2);
        }
        else
        {
            double t1 = b1->sig_t_vec[p1++];
            double t2 = b2->sig_t_vec[p2++];
            vector_len1 += (t1 * t1);
            vector_len2 += (t2 * t2);
            correlation += t1 * t2;
        }
    }
    while (p1 < b1->sig_t_vec.size())
    {
        long n1 = b1->sig_t_indx_vec[p1];
        double t1 = b1->sig_t_vec[p1++];
        vector_len1 += (t1 * t1);
    }
    while (p2 < b2->sig_t_vec.size())
    {
        long n2 = b2->sig_t_indx_vec[p2];
        double t2 = b2->sig_t_vec[p2++];
        vector_len2 += (t2 * t2);
    }

    return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

void CompareAllBacteria()
{
   Bacteria** b = new Bacteria*[number_bacteria];
   #pragma omp parallel for schedule(dynamic)
   for(int i=0; i<number_bacteria; i++)
	{
		printf("load %d of %d from %d\n", i+1, number_bacteria, omp_get_thread_num());
		b[i] = new Bacteria(bacteria_name[i]);
	}

//    #pragma omp parallel for
   for(int i=0; i<number_bacteria-1; i++)
      	#pragma omp parallel for schedule(dynamic, 2)
		for(int j=i+1; j<number_bacteria; j++)
		{
			printf("%2d %2d -> ", i, j);
			double correlation = CompareBacteria(b[i], b[j]);
			printf("%.20lf from %d\n", correlation, omp_get_thread_num());
		}
}

int main(int argc,char * argv[])
{
	auto start = std::chrono::high_resolution_clock::now();

	omp_set_num_threads(4);

	Init();
	ReadInputFile("../list.txt");
	CompareAllBacteria();

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;
	std::cout << "Comparision Time elapsed: " << elapsed.count() << " seconds\n";
	return 0;
}