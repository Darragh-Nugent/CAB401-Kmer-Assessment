#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <chrono>
#include <thread>
#include <vector>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>

int thread_count;

int number_bacteria;
char **bacteria_name;
long M_6, M_5, M_4;
short code[27] = {0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3};
#define encode(ch) code[ch - 'A']
#define LEN 6
#define AA_NUMBER 20
#define EPSILON 1e-010



using namespace std;

// Retrieved from https://www.geeksforgeeks.org/cpp/thread-pool-in-cpp/
class ThreadPool {
public:
    // // Constructor to creates a thread pool with given
    // number of threads
    ThreadPool(size_t num_threads
               = thread::hardware_concurrency())
    {

        // Creating worker threads
        for (size_t i = 0; i < num_threads; ++i) {
            threads_.emplace_back([this] {
                while (true) {
                    function<void()> task;
                    // The reason for putting the below code
                    // here is to unlock the queue before
                    // executing the task so that other
                    // threads can perform enqueue tasks
                    {
                        // Locking the queue so that data
                        // can be shared safely
                        unique_lock<mutex> lock(
                            queue_mutex_);

                        // Waiting until there is a task to
                        // execute or the pool is stopped
                        cv_.wait(lock, [this] {
                            return !tasks_.empty() || stop_;
                        });

                        // exit the thread in case the pool
                        // is stopped and there are no tasks
                        if (stop_ && tasks_.empty()) {
                            return;
                        }

                        // Get the next task from the queue
                        task = move(tasks_.front());
                        tasks_.pop();
                    }

                    task();
                }
            });
        }
    }

    // Destructor to stop the thread pool
    ~ThreadPool()
    {
        {
            // Lock the queue to update the stop flag safely
            unique_lock<mutex> lock(queue_mutex_);
            stop_ = true;
        }

        // Notify all threads
        cv_.notify_all();

        // Joining all worker threads to ensure they have
        // completed their tasks
        for (auto& thread : threads_) {
            thread.join();
        }
    }

    // Enqueue task for execution by the thread pool
    void enqueue(function<void()> task)
    {
        {
            unique_lock<std::mutex> lock(queue_mutex_);
            tasks_.emplace(move(task));
        }
        cv_.notify_one();
    }

private:
    // Vector to store worker threads
    vector<thread> threads_;

    // Queue of tasks
    queue<function<void()> > tasks_;

    // Mutex to synchronize access to shared data
    mutex queue_mutex_;

    // Condition variable to signal changes in the state of
    // the tasks queue
    condition_variable cv_;

    // Flag to indicate whether the thread pool should stop
    // or not
    bool stop_ = false;
};


void Init()
{
	M_4 = 1;
	for (int i = 0; i < LEN - 2; i++) // M_4 = AA_NUMBER ^ (LEN-2);
		M_4 *= AA_NUMBER;
	M_5 = M_4 * AA_NUMBER; // M_5 = AA_NUMBER ^ (LEN-1);
	M_6 = M_5 * AA_NUMBER; // M_6  = AA_NUMBER ^ (LEN);
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
	bacteria_name = new char *[number_bacteria];

	for (long i = 0; i < number_bacteria; i++)
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

void CreateBacteria(Bacteria **b, int thread_id, int start, int end)
{
	for (int i = start; i < end; i++)
	{
		printf("thread %d loading bacteria %d of %d\n", thread_id, i + 1, number_bacteria);
		b[i] = new Bacteria(bacteria_name[i]);
	}
}

void CompareBacteriaThread(Bacteria **b, int thread_id, int start, int end)
{
	for (int i = start; i < end; i++)
	{
		for (int j = i + 1; j < number_bacteria; j++)
		{
			printf("thread %d: %2d %2d -> ", thread_id, i, j);
			double correlation = CompareBacteria(b[i], b[j]);
			printf("%.20lf\n", correlation);
		}
	}
}

void CompareAllBacteria()
{
	Bacteria **b = new Bacteria *[number_bacteria];
	ThreadPool pool(thread_count);
	for (int i = 0; i < number_bacteria; i++)
	{
        pool.enqueue([b, i]() {
            printf("loading bacteria %d of %d\n", i + 1, number_bacteria);
            b[i] = new Bacteria(bacteria_name[i]);
        });
	}

	for (int i = 0; i < number_bacteria; i++)
	{
        pool.enqueue([b, i]() {
            for (int j = i + 1; j < number_bacteria; j++)
            {
                printf("comparing: %2d %2d -> ", i, j);
                double correlation = CompareBacteria(b[i], b[j]);
                printf("%.20lf\n", correlation);
            }
        });
	}

	delete[] b;
}

int main(int argc, char *argv[])
{
	time_t t1 = time(NULL);

	thread_count = 2;

	Init();
	ReadInputFile("../list.txt");
	// std::thread* threads;
	// CreateThreads(threads);
	CompareAllBacteria();

	time_t t2 = time(NULL);
	printf("time elapsed: %ld mers_5s\n", t2 - t1);
	return 0;
}