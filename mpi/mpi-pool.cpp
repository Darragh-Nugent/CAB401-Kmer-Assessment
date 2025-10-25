#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <chrono>
#include <iostream>

int rank, size;
const int tag_count = 1;
const int tag_vec = 2;
const int tag_vec_indx = 3;
const int tag_indx = 4;

int NAME_SIZE = 30;

int number_bacteria;
char *bacteria_name;
long M_6, M_5, M_4;
short code[27] = {0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};
#define encode(ch) code[ch - 'A']
#define LEN 6
#define AA_NUMBER 20
#define EPSILON 1e-010

void Init()
{
    M_4 = 1;
    for (int i = 0; i < LEN - 2; i++) // M_4 = AA_NUMBER ^ (LEN-2);
        M_4 *= AA_NUMBER;
    M_5 = M_4 * AA_NUMBER; // M_5 = AA_NUMBER ^ (LEN-1);
    M_6 = M_5 * AA_NUMBER; // M_6  = AA_NUMBER ^ (LEN);
}

struct BacteriaSummary
{
    int count;
    double *sig_t_vec;
    long *sig_t_indx_vec;
};

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
    double *sig_t_vec;
    long *sig_t_indx_vec;

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
        double *t = new double[M_6];

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
                t[i] = (mers_6[i] - stochastic) / stochastic;
                count++;
            }
            else
                t[i] = 0;
        }

        delete mers_5_div_total;
        delete mers_6;
        delete mers_5;

        sig_t_vec = new double[count];
        sig_t_indx_vec = new long[count];

        int pos = 0;
        for (long i = 0; i < M_6; i++)
        {
            if (t[i] != 0)
            {
                sig_t_vec[pos] = t[i];
                sig_t_indx_vec[pos] = i;
                pos++;
            }
        }
        delete t;

        fclose(bacteria_file);
    }
};

void ReadInputFile(const char *input_name)
{
    FILE *input_file = fopen(input_name, "r");
    if (input_file == NULL)
    {
        fprintf(stderr, "Error: failed to open file %s\n", input_name);
        exit(1);
    }

    fscanf(input_file, "%d", &number_bacteria);
    bacteria_name = new char[number_bacteria * NAME_SIZE];

    for (long i = 0; i < number_bacteria; i++)
    {
        char name[10];
        fscanf(input_file, "%s", name);
		snprintf(&bacteria_name[i * NAME_SIZE], NAME_SIZE, "../data/%s.faa", name);
    }
    fclose(input_file);
}

double CompareBacteria(BacteriaSummary *b1, BacteriaSummary *b2)
{
    double correlation = 0;
    double vector_len1 = 0;
    double vector_len2 = 0;
    long p1 = 0;
    long p2 = 0;
    while (p1 < b1->count && p2 < b2->count)
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
    while (p1 < b1->count)
    {
        long n1 = b1->sig_t_indx_vec[p1];
        double t1 = b1->sig_t_vec[p1++];
        vector_len1 += (t1 * t1);
    }
    while (p2 < b2->count)
    {
        long n2 = b2->sig_t_indx_vec[p2];
        double t2 = b2->sig_t_vec[p2++];
        vector_len2 += (t2 * t2);
    }

    return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

void RetrieveBacteriaInfo(int *all_counts, Bacteria **local_b, double **all_sig_t_vec, long **all_sig_t_indx_vec)
{
    for (int i = 0; i < number_bacteria; i++)
    {
        int owner_rank = i % size;

        if (rank == 0)
        {
            if (owner_rank == 0)
            {
                all_counts[i] = local_b[i]->count;
                all_sig_t_vec[i] = local_b[i]->sig_t_vec;
                all_sig_t_indx_vec[i] = local_b[i]->sig_t_indx_vec;
            }
            else
            {
                MPI_Recv(&all_counts[i], 1, MPI_INT, owner_rank, tag_count, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                all_sig_t_vec[i] = new double[all_counts[i]];
                MPI_Recv(all_sig_t_vec[i], all_counts[i], MPI_DOUBLE, owner_rank, tag_vec, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                all_sig_t_indx_vec[i] = new long[all_counts[i]];
                MPI_Recv(all_sig_t_indx_vec[i], all_counts[i], MPI_LONG, owner_rank, tag_vec_indx, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        else if (rank == owner_rank)
        {
            MPI_Send(&local_b[i]->count, 1, MPI_INT, 0, tag_count, MPI_COMM_WORLD);
            MPI_Send(local_b[i]->sig_t_vec, local_b[i]->count, MPI_DOUBLE, 0, tag_vec, MPI_COMM_WORLD);
            MPI_Send(local_b[i]->sig_t_indx_vec, local_b[i]->count, MPI_LONG, 0, tag_vec_indx, MPI_COMM_WORLD);
        }
    }
}

void CreateSummaries(int *all_counts, double **all_sig_t_vec, long **all_sig_t_indx_vec, BacteriaSummary **summaries)
{
    for (int i = 0; i < number_bacteria; i++)
    {
        MPI_Bcast(&all_counts[i], 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank != 0)
        {
            all_sig_t_vec[i] = new double[all_counts[i]];
            all_sig_t_indx_vec[i] = new long[all_counts[i]];
        }

        MPI_Bcast(all_sig_t_vec[i], all_counts[i], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(all_sig_t_indx_vec[i], all_counts[i], MPI_LONG, 0, MPI_COMM_WORLD);

        summaries[i] = new BacteriaSummary();
        summaries[i]->count = all_counts[i];
        summaries[i]->sig_t_vec = all_sig_t_vec[i];
        summaries[i]->sig_t_indx_vec = all_sig_t_indx_vec[i];
    }
}

void FreeBacteria(BacteriaSummary **summaries, int *all_counts, double **all_sig_t_vec, long **all_sig_t_indx_vec)
{
    for (int i = 0; i < number_bacteria; i++)
    {
        delete summaries[i];
        // delete local_b[i];
    }
    delete[] summaries;
    // delete[] local_b;
    delete[] all_counts;
    delete[] all_sig_t_vec;
    delete[] all_sig_t_indx_vec;
}

void RecvBacteria(int *all_counts, double **all_sig_t_vec, long **all_sig_t_indx_vec)
{
    int next_index = 0;
    int active_workers = size - 1;

    for (int i = 1; i < size && next_index < number_bacteria; i++)
    {
        MPI_Send(&next_index, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        next_index++;
    }

    while (active_workers > 0)
    {
        int index;
        MPI_Status status;

        MPI_Recv(&index, 1, MPI_INT, MPI_ANY_SOURCE, tag_indx, MPI_COMM_WORLD, &status);

        MPI_Recv(&all_counts[index], 1, MPI_INT, status.MPI_SOURCE, tag_count, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        all_sig_t_vec[index] = new double[all_counts[index]];
        MPI_Recv(all_sig_t_vec[index], all_counts[index], MPI_DOUBLE, status.MPI_SOURCE, tag_vec, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        all_sig_t_indx_vec[index] = new long[all_counts[index]];
        MPI_Recv(all_sig_t_indx_vec[index], all_counts[index], MPI_LONG, status.MPI_SOURCE, tag_vec_indx, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (next_index < number_bacteria)
        {
            MPI_Send(&next_index, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);

            next_index++;
        }
        else
        {
            int stop_signal = -1;
            MPI_Send(&stop_signal, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
            active_workers--;
        }
    }
}

void SendBacteria()
{
    int index;

    // Initially request the first index
    MPI_Recv(&index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    while (index != -1)
    {
        printf("Rank %d loading bacteria %d of %d\n", rank, index, number_bacteria);

        // Build the Bacteria object
        Bacteria *b = new Bacteria(&bacteria_name[index * NAME_SIZE]);

        // Send result to master
        MPI_Send(&index, 1, MPI_INT, 0, tag_indx, MPI_COMM_WORLD);
        MPI_Send(&b->count, 1, MPI_INT, 0, tag_count, MPI_COMM_WORLD);
        MPI_Send(b->sig_t_vec, b->count, MPI_DOUBLE, 0, tag_vec, MPI_COMM_WORLD);
        MPI_Send(b->sig_t_indx_vec, b->count, MPI_LONG, 0, tag_vec_indx, MPI_COMM_WORLD);

        delete b;

        // Ask for next index
        MPI_Recv(&index, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void CreateBacteria(int *all_counts, double **all_sig_t_vec, long **all_sig_t_indx_vec)
{
    if (rank == 0)
    {
        RecvBacteria(all_counts, all_sig_t_vec, all_sig_t_indx_vec);
    }
    else
    {
        SendBacteria();
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void CompareAllBacteria()
{
    int *all_counts = new int[number_bacteria];
    double **all_sig_t_vec = new double *[number_bacteria];
    long **all_sig_t_indx_vec = new long *[number_bacteria];

    auto start = std::chrono::high_resolution_clock::now();

    CreateBacteria(all_counts, all_sig_t_vec, all_sig_t_indx_vec);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    if (rank == 0) std::cout << "Bacteria creation time elapsed: " << elapsed.count() << " seconds\n";

    BacteriaSummary **summaries = new BacteriaSummary *[number_bacteria];

    // RetrieveBacteriaInfo(all_counts, all_sig_t_vec, all_sig_t_indx_vec);
    CreateSummaries(all_counts, all_sig_t_vec, all_sig_t_indx_vec, summaries);

    for (int i = rank; i < number_bacteria - 1; i += size)
    {
        for (int j = i + 1; j < number_bacteria; j++)
        {
            printf("Rank %d: %2d %2d -> ", rank, i, j);
            double correlation = CompareBacteria(summaries[i], summaries[j]);
            printf("%.20lf\n", correlation);
        }
    }

    FreeBacteria(summaries, all_counts, all_sig_t_vec, all_sig_t_indx_vec);
}

int main(int argc, char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Init();
    printf("MPI started with %d processes\n", size);
    if (rank == 0)
    {
	    ReadInputFile("../list.txt");
    }

    MPI_Bcast(&number_bacteria, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0)
        bacteria_name = new char[number_bacteria * NAME_SIZE];

    MPI_Bcast(bacteria_name, number_bacteria * NAME_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    printf("rank %d of %d\n", rank, size);

    CompareAllBacteria();

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time elapsed: " << elapsed.count() << " seconds\n";
    MPI_Finalize();
    return 0;
}