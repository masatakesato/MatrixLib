#include <stdio.h>
#include <math.h>

#include "GaussElimination.h"
#include "IndependentComponentAnalysis.h"









void ICA_Jutten91(Matrix *y, Matrix *W, Matrix *x){// ICA Jutten and Herault (1991)

	// f(y) = tanh(y), g(y) = y^3
	//        or
	// f(y) = y, g(y) = y^3

	unsigned int i,j,k;

	Matrix *dW = newMatrix(x->Row,x->Row);// ΔW
	Matrix_init_Random(W);// 最初にWの要素をランダムに決める

	for(int f=0; f<10; f++){

		for(i=0; i<dW->Row*dW->Column; i++) *(dW->Elements +i) =0;
		/**************** y = Wx ****************/
		Matrix_mult(y,W,x);

		/**************** ΔW ∝ -f(y)g(y)T ****************/
		for(i=0; i<dW->Row; i++){
			for(j=0; j<dW->Column; j++){
			
				for(k=0; k<y->Column; k++)
					*(dW->Elements + i*dW->Column + j) -= *(y->Elements + i*y->Column +k) * pow(*(y->Elements + j*y->Column +k),3);
			
				*(dW->Elements + i*dW->Column +j) *= 0.0001;

			}// end of j loop
		}// end of i loop

		/**************** W ← W + ΔW ****************/
		Matrix_add(W,W,dW);
		
		printf("%lf\n", Matrix_norm1(dW));
		Matrix_display(W);
		
		//if(Matrix_norm1(dW) <1.0e-3 && f<32) return;

	}// end of loop

}


void ICA_Fast(Matrix *a){





}