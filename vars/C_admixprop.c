
#include "admixprop.h"
#include <math.h>

void Q_update(const uint8_t* G, const double* zetabeta, const double* zetagamma, const double* xi, double* new_var, long N, long L, long K)
{
    uint8_t genotype;
    long n, l, k;
    double normbeta, normgamma;

    // loop over samples
    for (n=0; n<N; n++) {

        // loop over loci
        for (l=0; l<L; l++) {

            normbeta = 0.0;
            normgamma = 0.0;
            genotype = G[n*L+l];
            printf("\nQ update:\nN = %ld, L = %ld, g = %d\n", n, l, genotype);

            // missing data do not contribute
            if (genotype!=3) {

                // compute normalization
                for (k=0; k<K; k++) {
		  normbeta += zetabeta[l*K+k] * xi[n*K+k];
		  normgamma += zetagamma[l*K+k] * xi[n*K+k];
                }

                // loop over populations
                for (k=0; k<K; k++) {

                  // compute new estimate of variational parameters
                  // Though strangely expressed, this is actually
                  // E[Z_a] + E[Z_b]
		  new_var[n*K+k] += (((double) (2-genotype) * zetagamma[l*K+k] / normgamma) + \
                        ((double) genotype * zetabeta[l*K+k] / normbeta)) * xi[n*K + k];


                  int idx = l * K + k;
                  printf("z_a[%d, %d, %d] = %f\n", n, l, k,
                    zetabeta[idx] * xi[n * K + k] / normbeta);
                  printf("z_b[%d, %d, %d] = %f\n", n, l, k,
                    zetagamma[idx] * xi[n * K + k] / normgamma);

                } // k loop
            } // data missing check
        } // l loop
    } // n loop
}
