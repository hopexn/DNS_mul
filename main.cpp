#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iostream>
#include <cmath>
#include "mpi.h"
#include "matrix.h"
#include "dns.h"

using namespace std;

int main(int argc, char **argv) {
    double time_start, time_end;    //保存开始时间与结束时间
    MeshInfo meshInfo;              //保存拓扑信息
    FILE *fp;

    //创建并初始化三个矩阵， A、B为输入， C为b并行计算结果， D为串行计算结果（用于校验）
    Matrix *A = NULL;
    Matrix *B = NULL;
    Matrix *C = NULL;
    Matrix *D = NULL;

    //初始化
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &meshInfo.num_procs);             //获取节点数
    MPI_Comm_rank(MPI_COMM_WORLD, &meshInfo.my_rank);               //获取本节点ID
    meshInfo.dim_size = (int) pow(meshInfo.num_procs, 1 / 3.0);     //根据节点信息计算

//    //生成测试样例
//    for (int i = 0; i < 6; i++) {
//        for (int j = 0; j < 6; j++) {
//            cout << i * 6 + j << "\t";
//        }
//        cout << endl;
//    }

    //方阵大小
    if (meshInfo.my_rank == 0) {
        fp = fopen("yin.txt", "r");
        if (fp == NULL) {
            cout << "Cannot open file yin.txt" << endl;
            exit(1);
        }
        //读入矩阵大小
        fscanf(fp, "%d", &meshInfo.matrix_size);

        //创建并读入矩阵A
        A = CreateMatrix(meshInfo.matrix_size, meshInfo.matrix_size);
        for (int i = 0; i < meshInfo.matrix_size; i++) {
            for (int j = 0; j < meshInfo.matrix_size; j++) {
                fscanf(fp, "%d", &A->data[i][j]);
            }
        }

        //创建并读入矩阵B
        B = CreateMatrix(meshInfo.matrix_size, meshInfo.matrix_size);
        for (int i = 0; i < meshInfo.matrix_size; i++) {
            for (int j = 0; j < meshInfo.matrix_size; j++) {
                fscanf(fp, "%d", &B->data[i][j]);
            }
        }
    }

    MPI_Bcast(&meshInfo.matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);   //告诉网络上所有节点矩阵的大小

    create_topology(&meshInfo);                     //获取整个拓扑信息

    time_start = MPI_Wtime();                       //记录开始时间
    C = MultiplyMatrixesDNS(A, B, &meshInfo);       //进行DNS乘法
    time_end = MPI_Wtime();                         //记录结束时间

    MPI_Barrier(MPI_COMM_WORLD);

    //主节点输出结果
    if (meshInfo.my_rank == 0) {
        cout << endl << "Matrix A:" << endl;
        PrintMatrix(A);
        cout << endl << "Matrix B:" << endl;
        PrintMatrix(B);
        cout << endl << "Matrix C ( A x B ): " << endl;
        PrintMatrix(C);

        cout << endl << "Time cost: " << time_end - time_start << "s" << endl;

        //计算串行结果
        D = CreateMatrix(meshInfo.matrix_size, meshInfo.matrix_size);
        for (int i = 0; i < meshInfo.matrix_size; i++) {
            for (int j = 0; j < meshInfo.matrix_size; j++) {
                D->data[i][j] = 0;
                for (int k = 0; k < meshInfo.matrix_size; k++) {
                    D->data[i][j] += A->data[i][k] * B->data[k][j];
                }
            }
        }

        //PrintMatrix(D);  //输出串行结果过

        //将串行结果与并行结果作比较
        int flag = true;
        for (int i = 0; i < meshInfo.matrix_size; i++) {
            for (int j = 0; j < meshInfo.matrix_size; j++) {
                if (D->data[i][j] != C->data[i][j]) {
                    flag = false;
                    break;
                }
            }
            if (!flag) {
                break;
            }
        }

        //显示正确与否
        if (flag) {
            cout << "Correct" << endl;
        } else {
            cout << "Incorrect" << endl;
        }
        //释放资源
        FreeMatrix(A);
        FreeMatrix(B);
        FreeMatrix(C);
        FreeMatrix(D);
    }

    MPI_Finalize();

    return 0;

}

