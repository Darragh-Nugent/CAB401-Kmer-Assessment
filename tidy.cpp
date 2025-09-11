#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

int number_bacteria;
char** bacteria_name;
long M_6, M_5, M4;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};
#define encode(ch)		code[ch-'A']
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010

void Init()
{
	M4 = 1;
	for (int i=0; i<LEN-2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
		M4 *= AA_NUMBER; 
	M_5 = M4 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
	M_6  = M_5 *AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);
}

class Bacteria
{
private:
	long* mers_5;
	long mers_1[AA_NUMBER];
	long indexs;
	long total;
	long total_l;
	long complement;

	void InitVectors()
	{
		mers_6 = new long [M_6];
		mers_5 = new long [M_5];
		memset(mers_6, 0, M_6 * sizeof(long));
		memset(mers_5, 0, M_5 * sizeof(long));
		memset(mers_1, 0, AA_NUMBER * sizeof(long));
		total = 0;
		total_l = 0;
		complement = 0;
	}

	void init_buffer(char* buffer)
	{
		complement++;
		indexs = 0;
		for (int i=0; i<LEN-1; i++)
		{
			short enc = encode(buffer[i]);
			mers_1[enc]++;
			total_l++;
			indexs = indexs * AA_NUMBER + enc;
		}
		mers_5[indexs]++;
	}

	void cont_buffer(char ch)
	{
		short enc = encode(ch);
		mers_1[enc]++;
		total_l++;
		long index = indexs * AA_NUMBER + enc;
		mers_6[index]++;
		total++;
		indexs = (indexs % M4) * AA_NUMBER + enc;
		mers_5[indexs]++;
	}

public:
	long* mers_6;

	Bacteria(char* filename)
	{
		FILE* bacteria_file;
		errno_t OK = fopen_s(&bacteria_file, filename, "r");

		if (OK != 0)
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
				while (fgetc(bacteria_file) != '\n'); // skip rest of line

				char buffer[LEN-1];
				fread(buffer, sizeof(char), LEN-1, bacteria_file);
				init_buffer(buffer);
			}
			else if (ch != '\n')
				cont_buffer(ch);
		}
		fclose (bacteria_file);
	}

	~Bacteria()
	{
		delete mers_6;
		delete mers_5;
	}

	double stochastic_compute(long i)
	{
		double p1 = (double)mers_5[i / AA_NUMBER] / (total + complement);
		double p2 = (double) mers_1[i % AA_NUMBER] / total_l;
		double p3 = (double)mers_5[i % M_5] / (total + complement);
		double p4 = (double) mers_1[i / M_5] / total_l;
		return total * (p1*p2 + p3*p4) / 2;
	}
};

void ReadInputFile(const char* input_name)
{
	FILE* input_file;
	errno_t OK = fopen_s(&input_file, input_name, "r");

	if (OK != 0)
	{
		fprintf(stderr, "Error: failed to open file %s (Hint: check your working directory)\n", input_name);
		exit(1);
	}

    fscanf_s(input_file,"%d",&number_bacteria);
	bacteria_name = new char*[number_bacteria];

	for(long i=0;i<number_bacteria;i++)
	{
		char name[10];
		fscanf_s(input_file, "%s", name, 10);	
		bacteria_name[i] = new char[20];
		sprintf_s(bacteria_name[i], 20, "data/%s.faa", name);
	}
	fclose(input_file);
}

double CompareBacteria(Bacteria* b1, Bacteria* b2)
{
	double correlation = 0;
	double vector_len1=0;
	double vector_len2=0;

	for(long i=0; i<M_6; i++)
	{
		double stochastic1 = b1->stochastic_compute(i);
		double t1;
		if (stochastic1 > EPSILON) 
			t1 = (b1->mers_6[i] - stochastic1) / stochastic1;
		else 
			t1=0;
		vector_len1 += (t1 * t1);

		double stochastic2 = b2->stochastic_compute(i);
        double t2;
		if (stochastic2 > EPSILON) 
			t2 = (b2->mers_6[i] - stochastic2) / stochastic2;
		else 
			t2 = 0;
		vector_len2 += (t2 * t2);

		correlation = correlation + t1 * t2;
	}

	return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

void CompareAllBacteria()
{
    for(int i=0; i<number_bacteria-1; i++)
	{
		Bacteria* b1 = new Bacteria(bacteria_name[i]);

		for(int j=i+1; j<number_bacteria; j++)
		{
			Bacteria* b2 = new Bacteria(bacteria_name[j]);
			double correlation = CompareBacteria(b1, b2);
			printf("%03d %03d -> %.10lf\n", i, j, correlation);
			delete b2;
		}
		delete b1;
	}
}

int main(int argc,char * argv[])
{
	time_t t1 = time(NULL);

	Init();
	ReadInputFile("list.txt");
	CompareAllBacteria();

	time_t t2 = time(NULL);
	printf("time elapsed: %lld seconds\n", t2 - t1); 
	return 0;
}