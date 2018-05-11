#include <vector>
#include "lu_decomposition.hpp"


    void decompose(std::vector<float> &A, const int N) {

        for (int m = 0; m < N - 1; m++) {
//
			const int mref = m * N;
//
            for (int n = m + 1; n < N ; n++) {
//
				const int nref = n * N;			
//
				A[nref + m] = A[nref + m] / A[m + mref];
//
                for (int k = m + 1; k < N; k++) {
//
					A[k + nref] = A[k + nref] - A[m + nref] * A[k + mref];
//
				}
//
			}
//
		}
	}

    void solv(std::vector<float> &LU, std::vector<float> &B, const int N){

        float sum = 0.0;

		B[0] = B[0] / 1.0;

        for (int m = 1; m < N; m++) {

			const int mref = m * N;

			sum = 0.0;

            for (int n = 0; n < m ; n++) {

				sum = sum + LU[n + mref] * B[n];

			}

			B[m] = (B[m] - sum);

		}
//
//		Backward substitution - Ux=y
//
		B[N-1] = B[N-1] / LU[(N-1) * N + (N-1)];
//
        for (int m = N - 2; m > -1; m--) {
//
			sum = 0.0;
//
			const int mref = m * N;
//
            for (int n = m + 1; n < N; n++) {
//
				sum = sum + LU[mref + n] * B[n];
//
			}
//
			B[m] = (B[m] - sum) / LU[mref + m];
//
		}
	}
