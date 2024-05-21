#include <cmath>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <conio.h>
#include <malloc.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

using namespace std;

void secante();
void cuadraticaconfuncion();
void rk2b();
void rk3();
void rk4a();
void rk4b();
void linearecta();
void gaussjordan();
void montante();
void jacobi();
void gaussseidel();
void trapezoidal();
void simpson13();
void newton38();
void puntofijo();
void puntofijo();
void biseccion();
void euler();
void falsaposicion();
void cramer();
void raicespoli();
void cubica();
void cuadratica();
void eliminaciongaussiana();
void puntofijo();
void inversa();

void newtonraphson(){
	float n,error,e,x1,x,f,df;
	int a;
	cout<<"Ingrese numero a encontrar su raiz cubica"<<endl;
	cin>>n;
	cout<<"Indique cantidad maxima de iteraciones"<<endl;
	cin>>a;
	cout<<"Indique error de convergencia"<<endl;
	cin>>error;
	x1=n;
	f=(pow(x1,3))-n;
	df=3*x1*x1;
	x=x1-(f/df);
	e = abs(x-x1);
	for(int i = 0; i < a; i++){
		if(e>error){
			x1=x;
			f=(pow(x1,3))-n;
			df=3*x1*x1;
			x=x1-(f/df);
			
			e = abs(x-x1);
		}	
	}
	cout<<" x = "<<x<<endl;
}

void interpolacionlineal(){
	
    	std::string s1;
        int n, i, j;
        float x[100],  fx[100],  xint,  fxint,  s;
       cout<<"\n"<<s1;
       cout<<"\n\n  \4 Ingresa el numero de parejas de puntos: ";
       cin >> n;
       cout<<"\n"<<s;
       cout << "\n\n\nINGRESE VALORES DE x[i] y sus fx[i]\n\n";

		for ( i=0;  i<n;  i++ ) {

      cout << "\nx[ " << i << " ] : ";
      cin >> x[i];
      cout << "\nFx[ " << i << " ] : ";
      cin >> fx[i];
  }
      cout << "\n\nINGRESE VALOR A INTERPOLAR : " ;
      cin >> xint;

      fxint = 0;
	for( i=0; i<n; i++){
      s = 1;
      j = 0;
      while( j<n )
      {
             if( i != j){
                  s = s*( xint - x[j] )/( x[i] - x[j]);
                  }
                        j++;
            }
         fxint = fxint + s*fx[i];
    }
    cout<<"\n"<<s1;
cout << "\n\nEL VALOR DE LA FUNCION EN LA INTERPOLACION DE f( " << xint << ") es : " << fxint;
cout << endl << endl << endl;
}


void lagrange(){
	int n, x, y;
	float z,w,test,Yfinal,u,d;
	Yfinal=0;
	u=d=1;
	cout<<"Ingrese cantidad de pares a utilizar" << endl;
	cin>> n;
	float px[n],py[n];
	cout<<"Ingrese los valores de x:"<<endl;
	for (x = 0; x < n; x++){ //llena valores en x
			cout<< "x" <<x<<" :";
			cin>>px[x];
	}
	cout<<"Ingrese los valores de y:"<<endl;
	for (x = 0; x < n; x++){ //llena valores en y
			cout<< "y " <<x<<" :";
			cin>>py[x];
	}
	cout<<"Ingrese punto a evaluar:"<<endl;
	cin>>test;
	
	for (x = 0; x<n ; x++){
		for (y = 0; y<n ; y++){
			if (y!=x){
				u = u*(test-px[y]);
				d = d*(px[x]-px[y]);
			}
			else{
			}
		}
		
	Yfinal = Yfinal + (py[x]*(u/d));
	u=d=1;	
	}
	cout<<"Y ="<<Yfinal;
	cout<<endl;
}

void Newtonhaciaatras(){
	int n, x, y,m;
	double h,w,test,Yfinal,a,dy,d;
	Yfinal=dy=0;
	double u;
	u=a=1;
	cout<<"Ingrese cantidad de pares a utilizar" << endl;
	cin>> n;
	double px[n],py[n],c[n];
	for (x = 0; x < n-1; x++){ //hace 0 C
			c[x]=0;
	}
	cout<<"Ingrese los valores de x:"<<endl;
	for (x = 0; x < n; x++){ //llena valores en x
			cout<< "x" <<x<<" :";
			cin>>px[x];
	}
	cout<<"Ingrese los valores de y:"<<endl;
	for (x = 0; x < n; x++){ //llena valores en y
			cout<< "y " <<x<<" :";
			cin>>py[x];
	}
	cout<<"Ingrese punto a evaluar:"<<endl;
	cin>>test;
	h =abs(px[1]-px[0]);
	c[0]=py[n-1];
	for (x=1 ; x<n ; x++){
		m=1;
		dy=0;
		u=1;
		u = u*(x);
		for (y = 1 ; y<=x ; y++){
			dy =dy + ((u)*(pow(-1,y))*(py[n-y-1]));
			w=(x-y);
			d=(y+1);
			u = u*(w/d);
		}
		dy = dy + py[n-1];
		for(y=1;y<=x;y++){
			m=m*y;
		}
		c[x]=(dy)/(m*pow(h,x));
	}
	for (x = 0; x < n; x++){ 
			cout<< "c["<<x<<"] =" <<c[x]<<endl;
	}
	
	
	for (x=1; x<n;x++){
		for (y=1 ; y<=x ; y++){
			a = a*(test-px[n-y]);
		}
		a = a * c[x];
		Yfinal= Yfinal+a;
		a=1;
	}
	Yfinal = Yfinal + c[0];
	cout<<" Y = " <<Yfinal<<endl;
	
	
}
void Newtonhaciaadelante(){
	int n, x, y,m;
	double h,w,test,Yfinal,a,dy,d;
	Yfinal=dy=0;
	double u;
	u=a=1;
	cout<<"Ingrese cantidad de pares a utilizar" << endl;
	cin>> n;
	double px[n],py[n],c[n];
	for (x = 0; x < n-1; x++){ //hace 0 C
			c[x]=0;
	}
	cout<<"Ingrese los valores de x:"<<endl;
	for (x = 0; x < n; x++){ //llena valores en x
			cout<< "x" <<x<<" :";
			cin>>px[x];
	}
	cout<<"Ingrese los valores de y:"<<endl;
	for (x = 0; x < n; x++){ //llena valores en y
			cout<< "y " <<x<<" :";
			cin>>py[x];
	}
	cout<<"Ingrese punto a evaluar:"<<endl;
	cin>>test;
	h =abs(px[1]-px[0]);
	c[0]=py[0];
	for (x=1 ; x<n ; x++){
		m=1;
		dy=0;
		u=1;
		u = u*(x);
		for (y = 1 ; y<=x ; y++){
			dy =dy + ((u)*(pow(-1,y))*(py[x-y]));
			w=(x-y);
			d=(y+1);
			u = u*(w/d);
		}
		cout<<"dy = " <<dy<<endl;
		dy = dy + py[x];
		for(y=1;y<=x;y++){
			m=m*y;
		}
		
		c[x]=(dy)/(m*pow(h,x));
	}
	for (x = 0; x < n; x++){ 
			cout<< "c["<<x<<"] =" <<c[x]<<endl;
	}
	
	
	
	for (x=1; x<n;x++){
		for (y=0 ; y<x ; y++){
			a = a*(test-px[y]);
		}
		a = a * c[x];
		Yfinal= Yfinal+a;
		a=1;
	}
	Yfinal = Yfinal + c[0];
	cout<<" Y = " <<Yfinal<<endl;
}


int main() {
	cout <<"Ingrese que metodo quiere\n";
	cout <<"\n1 Newton Rhapson \n";
	cout <<"\n2 Lagrange\n";
	cout <<"\n3 Interpolacion lineal\n";
	cout <<"\n4 Gauss Jordan\n";
	cout <<"\n5 Montante\n";
	cout <<"\n6 Jacobi\n";
	cout <<"\n7 Gauss Seidel\n";
	cout<<"\nResolver la integral de (x^2 * e*x)dx por el metodo de :\n"<<endl;
	cout <<"\n8 Regla trapezoidal \n";
	cout <<"\n9 Simpson 1/3\n";
	cout <<"\n10 Simpson 3/8\n";
	cout <<"\n--------------------------------------------------------------\n";
	cout <<"\n11 Newton hacia atras\n";
	cout <<"\n12 Newton hacia adelante\n";
	cout <<"\n13 Biseccion\n";
	cout <<"\n14 Falsa posicion\n";
	cout <<"\n15 Linea recta\n";
	cout<<"\t\n  ---- Programa para resolver la ecuacion diferencial x*y ----"<<endl;
	cout <<"\n16 Metodo de Euler\n";
	cout <<"\n17 Runge-Kuta 2do orden\n";
	cout <<"\n18 Runge-Kuta Tercer orden\n";
	cout <<"\n19 Runge-Kuta 4to orden 1/3 de Simpson\n";
	cout <<"\n20 Runge-Kuta 4to orden 3/8 de Simpson\n";
	cout <<"\n--------------------------------------------------------------\n";
	cout <<"\n21 Cuadratica con funcion\n";
	cout <<"\n22 Cuadratica\n";
	cout <<"\n23 Raices de polinomio\n";
	cout <<"\n24 Cotes cerradas \n";
	cout <<"\n25 Cotes abiertas\n";
	cout <<"\n26 Cubica\n";
	cout <<"\n27 Lineal con funcion\n";
	cout <<"\n28 Eliminacion Gaussiana\n";
	cout <<"\n29 Punto fijo\n";
	
	
	int metodo = 0;
	cin >> metodo;
	switch(metodo){
		
		case 1: cout<<"\nNetwon Rhapson\n";
		newtonraphson();
		break;
		
		case 2: cout<<"\nLagrange\n";
		lagrange();
		break;
		
		case 3: cout<<"\nInterpolacion Lineal\n";
		interpolacionlineal();
		break;
		
		case 4: cout<<"\nGauss Jordan\n";
		gaussjordan();
		break;
		
		case 5: cout<<"\nMontante\n";
		montante();
		break;
		
		case 6: cout<<"\nJacobi\n";
		jacobi();
		break;
		
		case 7: cout<<"\nGauss Seidel\n";
		gaussseidel();
		break;
		
		case 8: cout<<"\nRegla trapezoidal\n";
		trapezoidal();
		break;
		
		case 9: cout<<"\nSimpson 1/3\n";
		simpson13();
		break;
		
		case 10: cout<<"\nSimpson 3/8\n";
		newton38();
		break;
		
		case 11: cout<<"\nNewton hacia atras\n";
		Newtonhaciaatras();
		break;
		
		case 12: cout<<"\nNewton hacia adelante\n";
		Newtonhaciaadelante();
		break;
		
		case 13: cout<<"\nBiseccion\n";
		biseccion();
		
		case 14: cout<<"\nFalsa Posicion\n";
		falsaposicion();
		break;
		
		case 15: cout<<"\nMetodo de Euler\n";
		euler();
		break;
		
		case 16: cout<<"\nMetodo de linea recta\n";
		linearecta();
		break;
		
		case 17: cout<<"\nRunge-Kuta 2do orden\n";
		rk2b();
		break;
		
		case 19: cout<<"\nRunge-Kuta Tercer orden\n";
		rk3();
		break;
		
		case 20: cout<<"\nRunge-Kuta 4to orden 1/3 de Simpson\n";
		rk4a();
		break;
		
		case 21: cout<<"\nRunge-Kuta 4to orden 3/8 de Simpson\n";
		rk4b();
		break;
		
		case 22: cout<<"\nCuadratica con funcion\n";
		cuadraticaconfuncion();
		break;
		
		case 23: cout<<"Cuadratica\n";
		cuadratica();
		break;
		
		
		case 24: cout<<"Raices de polinomios\n";
		raicespoli();
		break; 
		
		case 25: cout<<"Cotes cerradas\n";
		trapezoidal();
		break; 
		
		case 26: cout<<"Cotes abiertas\n";
		trapezoidal();
		break; 
		
		case 27: cout<<"Cubica\n";
		cubica();
		break;
		
		case 28: cout<<"Lineal con funcion\n";
		cuadraticaconfuncion();
		break;
		
		case 29: cout<<"Eliminacion Gaussiana\n";
		eliminaciongaussiana();
		break;
		
		case 30: cout<<"Punto fijo\n";
		puntofijo();
		break;
		
		
		default: cout<<"\nEstamos trabajando en ello\n";
		break;
	}
}



void cramer(){
	
	int n, x, y, i, j, a,p;
	p = 0;
	float z,w,q;
	float r=1;
	cout<<"Ingrese cantidad de variables (cantidad de ecuaciones)" << endl;
	cin>> n;
	
	float A[n][n];
	float AA[n][n];
	float B[n][n];
	float BB[n][n];
	float auxa[n];
	float mat1[n];
	float mat2[n];
	
		for (x = 0; x < n; x++){ // CREA 0 en aux
			auxa[x]=0;
			mat1[x]=0;
			mat2[x]=1;
		}
	for (x = 0; x < n; x++){ // lllena matriz identidad
		for(y = 0; y < n; y++){
			B[x][y]=0;
			BB[x][y]=0;
		}
		B[x][x]=1;
		BB[x][x]=1;
	}
	for (x = 0; x < n; x++){ //llena matriz A'
			cout<< "Ingrese los elementos del renglon " << x+1 <<" de la matriz de coeficientes"<<endl;
		for(y = 0; y < n; y++){
			cin>> A[x][y];
			AA[x][y]=A[x][y];
		}
	}
	cout<<"Ingrese la matriz de resultados"<<endl;
	for (x = 0; x<n ; x++){
		cin>>mat1[x];
	}
	
	cout<<"MATRIZ INGRESADA"<<endl;
	for (x = 0; x < n; x++){ //MUESTRA LA MATRIZ INGRESADA
		for(y = 0; y < n; y++){
			cout<<AA[x][y]<<"\t";
		}
		cout<<"| "<<mat1[x]<<endl;
	}
	
	for (x = 0; x < n; x++){ // Cambia renglones tal que no hayan 0 en la diagonal
		if (A[x][x] == 0){
			for (y = x; y < n ; y++){
				if (A[y][x] == 0){
				}
				else {
					for (a = 0; a < n; a++){
					auxa[a]=A[y][a];
					A[y][a]=A[x][a];
					A[x][a]=auxa[a];
				}
				B[x][x]=B[x][x] * (-1);
				}
			}	
		}
	}
		for (x = 0; x<n;x++){
		z=A[x][x];
		for (y = 0; y < n; y++){
			A[x][y]=A[x][y] / z; // CONSIGUE EL 1 EN LA DIAGONAL
		}
		B[x][x]=B[x][x]*z;
		for (i=x; i < n; i++){
			if(i == x){
			}
			else if (i != x){
				w = A[i][x];
				B[i][x]= w;	
				for (j = 0 ; j < n ; j++){	
					A[i][j]= (A[i][j])-(w * A[x][j]);
				
				}
			}
		
		}	
	}
		cout<<endl;
	for (x = 0; x < n ; x++){ // DETERMINANTE DE LA MATRIZ INGRESADA
		r =r*B[x][x];
	}
		for (x = 0; x < n ; x++){
			for(i = 0; i < n; i++){
				for(j = 0; j < n; j++){ //Cambia la columna de la matriz por la de los resultados
					A[i][j]=AA[i][j];
					B[i][j]=BB[i][j];
				}
				A[i][x]=mat1[i];
			}
			for (a = 0; a < n; a++){ // Cambia renglones tal que no hayan 0 en la diagonal
			if (A[a][a] == 0){
			for (y = a; y < n ; y++){
				if (A[y][a] == 0){
				}
				else {
					for (p = 0; p < n; p++){
					auxa[p]=A[y][p];
					A[y][p]=A[a][p];
					A[a][p]=auxa[p];
				}
				B[a][a]=B[a][a] * (-1);
				}
			}	
		}	
	}	
	for (a = 0; a<n;a++){
		z=A[a][a];
		for (y = 0; y < n; y++){
			A[a][y]=A[a][y] / z; // CONSIGUE EL 1 EN LA DIAGONAL
		}

		B[a][a]=B[a][a]*z;		
			for (i=a; i < n; i++){
			if(i == a){
			}
			else if (i != a){
				w = A[i][a];
				B[i][a]= w;	
				for (j = 0 ; j < n ; j++){	
					A[i][j]= (A[i][j])-(w * A[a][j]);
				
				}
			}
		
		}
	} 
		w = mat2[x];
			for (a = 0; a < n ; a++){
		w=w*B[a][a];
		}
		mat2[x]=w;
		}		
	cout<< "La DETERMINANTE de la matriz ingresada es = " << r;
		cout<<endl;
	
	cout<<"Los resultados para las incognitas son: "<<"\n\n";
	for (x = 0; x<n;x++){
		cout<<"X"<<x<<" = "<< mat2[x]<<"/"<<r<<" = " <<(mat2[x]/r) <<endl;
	}
	cout<<endl;	
}

void gaussjordan(){
		int n, x, y, i, j, a,m,p;
	float z,w;
	cout<<"Ingrese la cantida de variables a encontrar (cantidad de ecuaciones)" << endl;
	cin>> n;
	
	float A[n][n];
	float B[n];
	float auxa[n];
	float auxb=0;
	
		for (x = 0; x < n; x++){ // CREA 0 en aux
			auxa[x]=0;;
		}
	
	for (x = 0; x < n; x++){ //llena matriz A
			cout<< "Ingrese los elementos del renglon " << x+1 <<" de la matriz de coeficientes"<<endl;
		for(y = 0; y < n; y++){
			cin>> A[x][y];
		}
	}
		cout<<"Ingrese la matriz de resultados: "<<endl;
	for (x = 0; x < n; x++){ // llena matriz B de resultados
		cin>>B[x];
	}
		cout<<"MATRIZ INGRESADA"<<endl;
	for (x = 0; x < n; x++){ //MUESTRA LA MATRIZ INGRESADA
		for(y = 0; y < n; y++){
			cout<<A[x][y]<<"\t";
		}
		cout<<"| "<<B[x]<<endl;
		}
	for (x = 0; x < n; x++){ //QUITA 0 DE LA DIAGONAL CAMBIANDO RENGLONES
		if (A[x][x] == 0){
			for (y = x; y < n ; y++){
				if (A[y][x] == 0){
				}
				else {
					for (a = 0; a < n; a++){
					auxa[a]=A[y][a];
					A[y][a]=A[x][a];
					A[x][a]=auxa[a];
				}
				auxb=B[y];
				B[y]=B[x];
				B[x]=auxb;
				}
			}
			
		}
	}
		
		for (x = 0; x<n;x++){ 
		z=A[x][x];
		B[x] = B[x] / z;
		for (y = 0; y < n; y++){
			A[x][y]=A[x][y] / z; // CONSIGUE EL 1 EN LA DIAGONAL
		}
		for (i=0; i < n; i++){ //CREA LOS 0 EXCEPTO EN LA DIAGONAL
			if(i == x){
			}
			else if (i != x){
				w = A[i][x];
				for (j = 0 ; j < n ; j++){	
					A[i][j]= (A[i][j])-(w * A[x][j]);
				}
				B[i]= (B[i])-(w * B[x]);
				
			}
		}	
	}
	
	cout<<"MATRIZ GAUSSIANA: "<<endl;
	for (x = 0; x < n; x++){ //MUESTRA LA MATRIZ INGRESADA
		for(y = 0; y < n; y++){
			cout<<A[x][y]<<"\t";
		}
		cout<<"| "<<B[x]<<endl;
	}
	
	cout<<"Los valores de las incognitas son: "<< endl;
	for(x=0;x<n;x++){
		cout<<"X"<<x<<" = "<<B[x]<<endl;
	}
}


void montante(){
	int n, x, y,w;
	cout<<"Ingrese el grado de la matriz" << endl;
	cin>> n;
	
	float A[n][n]; //guarda la matriz anterior
	float B[n][n]; 
	float mata[n]; //guarda matriz de resultados anterior
	float matb[n];

	for (x = 0; x < n; x++){ //llena matriz A
			cout<< "Ingrese los elementos del renglon " << x+1 <<" de la matriz de coeficientes"<<endl;
		for(y = 0; y < n; y++){
			cin>> A[x][y];
			B[x][y]=0;
		}
		matb[x]=0;
	}
	cout<<"Ingrese la matriz de resultados: "<<endl;
	for (x=0;x<n;x++){
		cin>>mata[x];
		matb[x]=mata[x];
	}
	cout<<"MATRIZ A"<<endl;
	for (x = 0; x < n; x++){ //MUESTRA LA MATRIZ INGRESADA
		for(y = 0; y < n; y++){
			cout<<A[x][y]<<"\t";
		}
		cout<<endl;
	}
	for(x = 0; x<n ; x++){
		for (y = 0; y<n ; y++){ //iguala el renglon k de la siguiente matriz al de la pasada
			B[x][y]=A[x][y];
		}
		for (y = 0; y<n ; y++){ //crea 0 en la columna menos en la diagonal
			if (y == x){
			}
			else {
				B[y][x]=0;
			}
		}
		for (y = 0; y<n;y++){ //iguala el valor pasado de la diagonal
			B[y][y]=A[x][x];
		}
		for (y = 0; y<n; y++){ //hace el producto cruz
			if (y == x){
			}
			else {
				matb[y]=(mata[y]*A[x][x]) - (A[y][x]*mata[x]);
				if (x == 0){
				}
				else {
					matb[y]=matb[y] / A[0][0];
				}
			}
			for (w = x+1 ; w<n ; w++){
			if (y == x){
			}
			else {
				B[y][w]=(A[y][w]*A[x][x]) - (A[y][x]*A[x][w]);
				if (x == 0){
				}
				else {
					B[y][w]=B[y][w] / A[0][0];
				}
			}
		}
		}
		for (y = 0; y < n ; y++){
			for (w = 0; w<n; w++){
				A[y][w]=B[y][w];																	
			}
			mata[y]=matb[y];
		}
	}	
	cout<<"Matriz final: "<<endl;
	for (x = 0; x < n; x++){ //MUESTRA LA MATRIZ INGRESADA
		for(y = 0; y < n; y++){
			cout<<B[x][y]<<"\t";
		}
		cout<<"| "<<matb[x]<<endl;
		cout<<endl<<endl;
	}
cout<<"Los valores de las variables son: \n\n";
for(x = 0; x<n ; x++){
	cout<<"X"<<x<<" = "<<matb[x] / B[x][x] <<endl;
}	
}


void jacobi(){
		int n, x, y,a,i,w;
	float error=.1;
	float z;
	i=1;
	cout<<"Para este programa es necesario que la matriz ingresada tenga los coeficientes mas altos en la diagonal"<<endl;
	cout<<"Ingrese la cantida de variables a encontrar (cantidad de ecuaciones)" << endl;
	cin>> n;
	
	float A[n][n];
	float B[n];
	float sig[n];
	float ant[n];
	float auxa[n];
	float auxb=0;
	
		for (x = 0; x < n; x++){ // CREA 0 en aux
			auxa[x]=0;
			sig[x]=0;
		}
	
	for (x = 0; x < n; x++){ //llena matriz A
			cout<< "Ingrese los elementos del renglon " << x+1 <<" de la matriz de coeficientes"<<endl;
		for(y = 0; y < n; y++){
			cin>> A[x][y];
		}
	}
		cout<<"Ingrese la matriz de resultados: "<<endl;
	for (x = 0; x < n; x++){ // llena matriz B de resultados
		cin>>B[x];
	}
	cout<<"Ingrese cantidad maxima de iteraciones"<<endl;
	cin>>w;
	cout<<"Ingrese el error de convergencia"<<endl;
	cin>>error;
		cout<<"MATRIZ INGRESADA"<<endl;
	for (x = 0; x < n; x++){ //MUESTRA LA MATRIZ INGRESADA
		for(y = 0; y < n; y++){
			cout<<A[x][y]<<"\t";
		}
		cout<<"| "<<B[x]<<endl;
		}
	for (x = 0; x < n; x++){ //QUITA 0 DE LA DIAGONAL CAMBIANDO RENGLONES
		if (A[x][x] == 0){
			for (y = x; y < n ; y++){
				if (A[y][x] == 0){
				}
				else {
					for (a = 0; a < n; a++){
					auxa[a]=A[y][a];
					A[y][a]=A[x][a];
					A[x][a]=auxa[a];
				}
				auxb=B[y];
				B[y]=B[x];
				B[x]=auxb;
				}
			}
			
		}
	}
	for(x = 0; x<n ; x++){
		z=A[x][x];
		for(y = 0; y<n ; y++){
		A[x][y]=(A[x][y])/(-z);
		}
		A[x][x]=0;
		B[x]=B[x] / z;
		ant[x]=B[x];
	}
	for(a = 0; a<n ; a++){
		for (y = 0; y<n ; y++){
		sig[a]=sig[a]+(ant[y]*A[a][y]);
		}
		sig[a]=sig[a]+B[a];
		}
	for (x = 0; x<n ; x++){
		while (abs(sig[x]-ant[x]) > error && w > i){ //Hace los pasos hasta que el error sea menor a error
			for (y = 0; y<n ;y++){
				ant[y]=sig[y];
			}
		for(a = 0; a<n ; a++){
				sig[a]=0;
				for (y = 0; y<n ; y++){
				sig[a]=sig[a]+(ant[y]*A[a][y]);
			}
			sig[a]=sig[a]+B[a];
		}
		i++;	
		}
}
if(i > w ){
	cout<<"El error no converge a "<<error<<" en "<<w<<" iteraciones"<<endl;
}
else{
cout<<endl<<" Se llego a los resultados en "<<i+1<<" iteraciones\n\n";
cout<<" Con un error de "<<error<<endl;
cout<<" Los valores de las variables son: \n\n";
for(x = 0; x<n ; x++){
	cout<<"X"<<x<<" = "<<sig[x]<</*"	Converge a = "<<round(sig[x])<<*/endl;
}
}
}


void gaussseidel(){
	int n, x, y,a,i,w;
	float error=.01;
	float z;
	i = 1;
	cout<<"Para este programa es necesario que la matriz ingresada tenga los coeficientes mas altos en la diagonal"<<endl;
	cout<<"Ingrese la cantida de variables a encontrar (cantidad de ecuaciones)" << endl;
	cin>> n;
	
	float A[n][n];
	float B[n];	
	float sig[n];
	float ant[n];
	float auxa[n];
	float auxb=0;
	float b[n];
	
		for (x = 0; x < n; x++){ // CREA 0 en aux
			auxa[x]=0;
			sig[x]=0;
			ant[x]=0;
			b[x]=0;
		}
	
	for (x = 0; x < n; x++){ //llena matriz A
			cout<< "Ingrese los elementos del renglon " << x+1 <<" de la matriz de coeficientes"<<endl;
		for(y = 0; y < n; y++){
			cin>> A[x][y];
		}
	}
		cout<<"Ingrese la matriz de resultados: "<<endl;
	for (x = 0; x < n; x++){ // llena matriz B de resultados
		cin>>B[x];
	}
	cout<<"Ingrese cantidad maxima de iteraciones"<<endl;
	cin>>w;
	cout<<"Ingrese error de convergencia"<<endl;
	cin>>error;
		cout<<"MATRIZ INGRESADA"<<endl;
	for (x = 0; x < n; x++){ //MUESTRA LA MATRIZ INGRESADA
		for(y = 0; y < n; y++){
			cout<<A[x][y]<<"\t";
		}
		cout<<"| "<<B[x]<<endl;
		}
	for (x = 0; x < n; x++){ //QUITA 0 DE LA DIAGONAL CAMBIANDO RENGLONES
		if (A[x][x] == 0){
			for (y = x; y < n ; y++){
				if (A[y][x] == 0){
				}
				else {
					for (a = 0; a < n; a++){
					auxa[a]=A[y][a];
					A[y][a]=A[x][a];
					A[x][a]=auxa[a];
				}
				auxb=B[y];
				B[y]=B[x];
				B[x]=auxb;
				}
			}
			
		}
	}

	for(x = 0; x<n ; x++){
		z=A[x][x];
		for(y = 0; y<n ; y++){
		A[x][y]=(A[x][y])/(-z);
		}
		A[x][x]=0;
		B[x]=B[x] / z;
	}
	for(a = 0; a<n ; a++){ //hace el primer paso
		for (y = 0; y<n ; y++){
		sig[a]=sig[a]+(ant[y]*A[a][y]);
		}
		sig[a]=sig[a]+B[a];
		ant[a]=sig[a];
		}
	for (x = 0; x<n ; x++){
		while (abs(sig[x]-b[x]) > error && w>i){ //Hace los pasos hasta que el error sea menor a error
			for (y = 0; y<n ;y++){
				b[y]=sig[y];
			}
		for(a = 0; a<n ; a++){
				sig[a]=0;
				for (y = 0; y<n ; y++){
				sig[a]=sig[a]+(ant[y]*A[a][y]);
			}
			sig[a]=sig[a]+B[a];
			ant[a]=sig[a];
		}
		i++;	
		}
}
if(i > w){
	cout<<"No el error no converge a "<<error<<" en "<<w<<" iteraciones"<<endl;
}
else{
cout<<endl<<" Se llego a los resultados en "<<i+1<<" iteraciones\n\n";
cout<<" Los valores de las variables son: \n\n";
for(x = 0; x<n ; x++){
	cout<<"X"<<x<<" = "<<sig[x]<</*"	Converge a = "<<round(sig[x])<<*/endl;
}
}
}

void trapezoidal(){
	float xi,xf,h,w,i;
w=0;
cout<<"Ingrese limite inferior de la integral : "<<endl;
cin>>xi;
cout<<"Ingrese limite superior de la integral : "<<endl;
cin>>xf;
cout<<"h = ";
cin>>h;

for (i=xi;i<=xf;i=i+h){
	if (i ==xi || i == xf ){
		w = w + (pow(i,2)*exp(i)*(h/2));
	}
	else{
		w = w + (pow(i,2)*exp(i)*(h));
	}
}
cout<<"El resultado obtenido es = "<<w<<endl;
cout<<"El resultado exacto es = "<<(exp(xf)*(pow(xf,2)-(2*xf)+2))-(exp(xi)*(pow(xi,2)-(2*xi)+2))<<endl;
cout<<"El error es = "<<abs(w-((exp(xf)*(pow(xf,2)-(2*xf)+2))-(exp(xi)*(pow(xi,2)-(2*xi)+2))))<<endl;
}
void simpson13(){
	float xi,xf,h,w,i;
int a=0;
w=0;
cout<<"Ingrese limite inferior de la integral : "<<endl;
cin>>xi;
cout<<"Ingrese limite superior de la integral : "<<endl;
cin>>xf;
cout<<"h = ";
cin>>h;

for (i=xi;i<=xf;i=i+h){
	if (i ==xi || i == xf ){
		w = w + (pow(i,2)*exp(i)*(h/3));
	}
	else{
		if (a%2 == 0){
		w = w + (pow(i,2)*exp(i)*(2*h/3));
	}
		else {
		w = w + (pow(i,2)*exp(i)*(4*h/3));	
		}
	}
	a++;
}

cout<<"El resultado obtenido es = "<<w<<endl;
cout<<"El resultado exacto es = "<<(exp(xf)*(pow(xf,2)-(2*xf)+2))-(exp(xi)*(pow(xi,2)-(2*xi)+2))<<endl;
cout<<"El error es = "<<abs(w-((exp(xf)*(pow(xf,2)-(2*xf)+2))-(exp(xi)*(pow(xi,2)-(2*xi)+2))))<<endl;
}
void newton38(){
	float xi,xf,h,w,i;
int a=0;
w=0;
cout<<"Ingrese limite inferior de la integral : "<<endl;
cin>>xi;
cout<<"Ingrese limite superior de la integral : "<<endl;
cin>>xf;
cout<<"h = ";
cin>>h;
for (i=xi;i<=xf;i=i+h){
	if (i ==xi || i+h >= xf ){
		w = w + (pow(i,2)*exp(i)*(3*h/8));
	}
	else{
		if (a%3 == 0){
		w = w + (pow(i,2)*exp(i)*(6*h/8));
	}
		else {
		w = w + (pow(i,2)*exp(i)*(9*h/8));	
		}
	}
	a++;
}
cout<<"El resultado obtenido es = "<<w<<endl;
cout<<"El resultado exacto es = "<<(exp(xf)*(pow(xf,2)-(2*xf)+2))-(exp(xi)*(pow(xi,2)-(2*xi)+2))<<endl;
cout<<"El error es = "<<abs(w-((exp(xf)*(pow(xf,2)-(2*xf)+2))-(exp(xi)*(pow(xi,2)-(2*xi)+2))))<<endl;
}


#include <iostream>
#include <iomanip> // setprecision
#include <cmath>

#define PRECISION 6

using namespace std;

double f(double x);
void imprimePuntos(double a, double b);

void biseccion()
{
   cout << setprecision(PRECISION); // Establecemos la precisión
   
   double a, b, tolerancia;
   
   cout << "\nCalculo de las raices de una funcion aplicando el metodo de la biseccion" << endl;
   cout << "\nIngrese el intervalo inicial [a, b]" << endl;
   cout << "\na = ";
   cin >> a;
   
   cout << "b = ";
   cin >> b;
   
   imprimePuntos(a, b);
   
   cout << "\nEscoja el intervalo adecuado" << endl;
   cout << "\na = ";
   cin >> a;
   
   cout << "b = ";
   cin >> b;
   
   // [a, b]
   float xr; // raiz de la función
   
   /*
   cout << "\nf(" << a << ") = " << f(a) << endl;
   cout << "f(" << b << ") = " << f(b) << endl;
   */
   
   if (f(a) * f(b) > 0) {
      
      cout << "\nNo se puede aplicar el metodo de la biseccion\n";
      cout << "porque f(" << a << ") y f(" << b << ") tienes el mismo signo" << endl;
      
   } else {
      cout << "Tolerancia = ";
      cin >> tolerancia;
      cout << "\na\tb\tx\tf(a)\t\tf(b)\t\tf(x)\n" << endl;
      do {
         xr = (a + b) / 2.0;
         cout << a << "\t" << b << "\t" << xr << "\t";
         cout << f(a) << "\t" << f(b) << "\t" << f(xr) << endl;

         // Vemos si cumple o no cumple
         if (abs(f(xr)) <= tolerancia) { // xr sería la raiz de f
         
            cout << "\n\nPara una tolerancia " << tolerancia << " la raiz de f es " << xr << endl;
         
            break;
            
         } else {
            if (f(xr) * f(a) > 0) {
               a = xr;
            } else if (f(xr) * f(b) > 0) {
               b = xr;
            }
         }
         
      } while (1);
   }
   
   cin.get();
   cin.get();
}


double f(double x) 
{
   return exp(-1 * x) - cos(3 * x) - 0.5;
}

#define INTERVALOS 10
void imprimePuntos(double a, double b)
{
   int puntos = INTERVALOS + 1;
   
   double ancho = (b - a) / INTERVALOS;
   
   cout << "\n";
   cout << "\tx\tf(x)\n" << endl;
   for (int i = 0; i < puntos; i++) {
      cout << "\t" << a << "\t" << f(a) << endl;
      a = a + ancho;
   }
}

#define PRECISION 6
#define INTERVALOS 10

using namespace std;

void tabula(double a, double b);
double f(double x);

void falsaposicion()
{
	cout << setprecision(PRECISION);
	cout << "\nCalculo de las raices de una funcion aplicando el metodo de la falsa posicion\n";
	cout << "\nIngrese el intervalo inicial [a,b]:" << endl;
	
	double a, b, tolerancia;
	
	cout << "\na = ";
	cin >> a;
	
	cout << "b = ";
	cin >> b;
	
	tabula(a, b);
	
	cout << "\nEscoja el intervalo adecuado" << endl;
	cout << "\na = ";
	cin >> a;
	
	cout << "b = ";
	cin >> b;
	
	double xr; // La solución aproximada
	double xa = 0; // Solución anterior
	
	if (f(a) * f(b) > 0) {
		cout << "\nNo se puede aplicar el metodo de la falsa posicion\n";
		cout << "porque f(" << a << ") y f(" << b << ") tienen el mismo signo" << endl;
	
	} else {
		cout << "Tolerancia = ";
		cin >> tolerancia;
		
		cout << "\na\tb\tx\tf(a)\t\tf(b)\t\tf(x)\n" << endl;
		do {
			xr = b - f(b) * ((b - a) / (f(b) - f(a)));
			
			cout << a << "\t" << b << "\t" << xr << "\t" << f(a) << "\t" << f(b) << "\t" << f(xr) << endl;
			
			// if (fabs(f(xr)) <= tolerancia) {
			if (fabs(xr - xa) / fabs(xr) <= tolerancia) {
				cout << "\n\nPara una tolerancia de " << tolerancia << " la raiz de f es: " << xr << endl;
				break;
			
			} else {
				xa = xr; // Se guarda el valor de la aproximación anterior
				if (f(xr) * f(a) > 0) {
					a = xr;
			
				} else if (f(xr) * f(b) > 0) {
					b = xr;
				}
			}
		
		} while (1);
	}
	
	cin.get();
	cin.get();
}

void tabula(double a, double b)
{
	int puntos = INTERVALOS + 1;
	
	double ancho = (b - a) / INTERVALOS;
	
	cout << "\n\tx\tf(x) " << endl;
	for (int i = 0; i < puntos; i++) {
		cout << "\t" << a << "\t" << f(a) << endl;
		a = a + ancho;
	}
}

float f(float x)
{
	// return exp(-1 * x) - cos(3 * x) - 0.5;
	return exp(-1 * x) - cos(x);
}

void euler(){
	int i,a;
	float xi,xf,h,k,y,f;
	
	cout<<"Ingrese valor inicial de x : \nx0 = ";
	cin>>xi;
	cout<<"Ingrese valor final de x : \nxf =";
	cin>>xf;
	cout<<"Ingrese valor inicial de y : \ny =";
	cin>>y;
	cout<<"Ingrese numero de iteraciones : \ni =";
	cin>>i;
	h = (xf - xi)/i;
	cout<<"\nh = "<<h<<endl;
	cout<<"0\t|x = "<<xi<<"\ty = "<<y;
	f = (xi*y);
	cout<<"\tf(x,y) = "<<f<<endl;
	k = xi+h;
	for(a = 1; a<=i ; a++)
	{
		y = y + (h*f);
		f = (k*y);
		cout<<a<<"\t|x = "<<k<<"\ty = "<<y<<"\tf(x,y) = "<<f<<endl;
		k=k+h;
	}
	cout<<"\nEl resultado es : \n\ty = "<<y<<endl;
}
  

void linearecta()
 {  double a,b,c,d,m;
    cout<<"Digite dos puntos P=(a,b) y Q(c,d):"<<endl<<endl;
    cout<<"a : "; cin>>a; cout<<"b : "; cin>>b;
    cout<<"c : "; cin>>c; cout<<"d : "; cin>>d;
    cout<<endl;
    if (b!=d) 
    {     m = (a-c)/(b-d);
          if (m>0)
          {cout<<"La recta determinada por estos dos puntos P=("<<a<<","<<b<<") y Q("<<c<<","<<d<<") es Creciente"<<endl;}
          else 
          {cout<<"La recta determinada por estos dos puntos P=("<<a<<","<<b<<") y Q("<<c<<","<<d<<") es Decreciente"<<endl;}
          if (m==0)
          {cout<<"La recta determinada por estos dos puntos P=("<<a<<","<<b<<") y Q("<<c<<","<<d<<") es Horizontal"<<endl;}
    } 
    else 
    {   if(a!=c)
        {cout<<"La recta determinada por estos dos puntos P=("<<a<<","<<b<<") y Q("<<c<<","<<d<<") es Vertical"<<endl;}
        else 
        {cout <<"Lo sentimos!. Por un mismo punto pasan infinitas rectas..."<<endl;}
    }
    cin.get(); /*Recuerda que esta linea es por si usas Windows*/
    cin.get();
 }



void rk2b(){
	int i,a;
	float xi,xf,h,k,y,k1,k2;
	
	cout<<"Ingrese valor inicial de x : \nx0 = ";
	cin>>xi;
	cout<<"Ingrese valor final de x : \nxf =";
	cin>>xf;
	cout<<"Ingrese valor inicial de y : \ny =";
	cin>>y;
	cout<<"Ingrese numero de iteraciones : \ni =";
	cin>>i;
	h = (xf - xi)/i;
	cout<<"\nh = "<<h<<endl;
	cout<<"0\t|x = "<<xi<<"\ty = "<<y;
	k1 = (xi*y);
	k2 = (xi+(h/2))*(y+((h/2)*k1));
	cout<<"\tk1 = "<<k1<<"\tk2 = "<<k2<<endl;
	k = xi+h;
	for(a = 1; a<=i ; a++)
	{
		y = y + (h*k2);
		k1 = (k*y);
		k2 = (k+(h/2))*(y+((h/2)*k1));
		cout<<a<<"\t|x = "<<k<<"\ty = "<<y<<"\tk1 = "<<k1<<"\tk2 = "<<k2<<endl;
		k=k+h;
	}
	cout<<"\nEl resultado es : \n\ty = "<<y<<endl;
}
void rk3(){
	int i,a;
	float xi,xf,h,k,y,k1,k2,k3;
	
	cout<<"Ingrese valor inicial de x : \nx0 = ";
	cin>>xi;
	cout<<"Ingrese valor final de x : \nxf =";
	cin>>xf;
	cout<<"Ingrese valor inicial de y : \ny =";
	cin>>y;
	cout<<"Ingrese numero de iteraciones : \ni =";
	cin>>i;
	h = (xf - xi)/i;
	cout<<"\nh = "<<h<<endl;
	cout<<"0   |x = "<<xi<<"   y = "<<y;
	k1 = (xi*y);
	k2 = (xi+(h/2))*(y+((h/2)*k1));
	k3 = (xi+h)*(y+(2*h*k2)-(h*k1));
	cout<<"   k1 = "<<k1<<"   k2 = "<<k2<<"   k3 = "<<k3<<endl;
	k = xi+h;
	for(a = 1; a<=i ; a++)
	{
		y = y + ((h/6)*(k1+(4*k2)+k3));
		k1 = (k*y);
		k2 = (k+(h/2))*(y+((h/2)*k1));
		k3 = (k+h)*(y+(2*h*k2)-(h*k1));
		cout<<a<<"   |x = "<<k<<"   y = "<<y<<"   k1 = "<<k1<<"   k2 = "<<k2<<"   k3 = "<<k3<<endl;
		k=k+h;
	}
	cout<<"\nEl resultado es : \n\ty = "<<y<<endl;
}
void rk4a(){
	int i,a;
	float xi,xf,h,k,y,k1,k2,k3,k4;
	
	cout<<"Ingrese valor inicial de x : \nx0 = ";
	cin>>xi;
	cout<<"Ingrese valor final de x : \nxf =";
	cin>>xf;
	cout<<"Ingrese valor inicial de y : \ny =";
	cin>>y;
	cout<<"Ingrese numero de iteraciones : \ni =";
	cin>>i;
	h = (xf - xi)/i;
	cout<<"\nh = "<<h<<endl;
	cout<<"0  |x= "<<xi<<"  y= "<<y;
	k1 = (xi*y);
	k2 = (xi+(h/2))*(y+((h/2)*k1));
	k3 = (xi+(h/2))*(y+((h/2)*k2));
	k4 = (xi+h)*(y+(h*k3));
	
	cout<<"  k1= "<<k1<<"  k2= "<<k2<<"  k3= "<<k3<<"  k4= "<<k4<<endl;
	k = xi+h;
	for(a = 1; a<=i ; a++)
	{
		y = y + ((h/6)*(k1+(2*k2)+(2*k3)+k4));
		k1 = (k*y);
		k2 = (k+(h/2))*(y+((h/2)*k1));
		k3 = (k+(h/2))*(y+((h/2)*k2));
		k4 = (k+h)*(y+(h*k3));
		
		cout<<a<<"  |x= "<<k<<"  y= "<<y<<"  k1= "<<k1<<"  k2= "<<k2<<"  k3= "<<k3<<"  k4= "<<k4<<endl;
		k=k+h;
	}
	cout<<"\nEl resultado es : \n\ty = "<<y<<endl;
}
void rk4b(){
	
	int i,a;
	float xi,xf,h,k,y,k1,k2,k3,k4;
	
	cout<<"Ingrese valor inicial de x : \nx0 = ";
	cin>>xi;
	cout<<"Ingrese valor final de x : \nxf =";
	cin>>xf;
	cout<<"Ingrese valor inicial de y : \ny =";
	cin>>y;
	cout<<"Ingrese numero de iteraciones : \ni =";
	cin>>i;
	h = (xf - xi)/i;
	cout<<"\nh = "<<h<<endl;
	cout<<"0  |x= "<<xi<<"  y= "<<y;
	k1 = (xi*y);
	k2 = (xi+(h/3))*(y+((h/3)*k1));
	k3 = (xi+(2*h/3))*(y+(h*k2)-((h/3)*k1));
	k4 = (xi+h)*(y+(h*k3)-(h*k2)+(h*k1));
	
	cout<<"  k1= "<<k1<<"  k2= "<<k2<<"  k3= "<<k3<<"  k4= "<<k4<<endl;
	k = xi+h;
	for(a = 1; a<=i ; a++)
	{
		y = y + ((h/8)*(k1+(3*k2)+(3*k3)+k4));
		k1 = (k*y);
		k2 = (k+(h/3))*(y+((h/3)*k1));
		k3 = (k+(2*h/3))*(y+(h*k2)-((h/3)*k1));
		k4 = (k+h)*(y+(h*k3)-(h*k2)+(h*k1));
		
		cout<<a<<"  |x= "<<k<<"  y= "<<y<<"  k1= "<<k1<<"  k2= "<<k2<<"  k3= "<<k3<<"  k4= "<<k4<<endl;
		k=k+h;
	}
	cout<<"\nEl resultado es : \n\ty = "<<y<<endl;
}


void cuadraticaconfuncion() 
{
float a,b,c,d,x,y;
cout << "ingrese a:";
cin >> a;
cout << "ingrese b:";
cin >> b;
cout << "ingrese c:";
cin >> c;
if (a == 0)
cout << "no cumple la condición."<<"\n";
else
{
d = (b*b) - (4*a*c);
if (d == 0)
{
x = (-b/2*a);
y = x;
}
}
	cout << "y="<<y;
}

void cuadratica()

{

float a,b,c,d,x1,x2;

printf ("Listo para encontrar las soluciones de una ecuacion cuadratica?");

printf ("\t\t\ttIntroduce el coeficiente del termino cuadratico:");

scanf ("%f",&a);

printf ("\tIntroduce el coeficiente del termino lineal:");

scanf ("%f",&b);

printf ("\tIntroduce el coeficiente del termino independiente:");

scanf ("%f",&c);

if (a!=0){

printf ("Aqui vamos...");

}else{


printf ("No es posible realizar la operacion"); }

{

d=sqrt(b*b-(4*a*b));

}

if (d>0)

{

x1=((b*-1)+(d))/(2*a);

x2=((b*-1)-(d))/(2*a);

printf ("\t\El resultado de x1 es: %f",x1);

printf ("\t\El resultado de x2 es: %f",x2);

}else{


printf("No es posible realizar la operacion, revisa tus datos");

getch ();

}
}

void raicespoli()
{

	float  a, b ,c, rta1, rta2, lin, base, raiz, raiz2, img, x;
    cout<<"\n                                 ax^2 + bx + c = 0"<<endl;
	cout<<endl;
	cout<<"Ingrese el valor de a:";
	cin>>a;
	cout<<endl;
	cout<<"Ingrese el valor de b:";
	cin>>b;
	cout<<endl;
	cout<<"Ingrese el valor de c:";
	cin>>c;
	cout<<endl;
	system("cls");
	cout<<endl;
        cout<<a<<"x^2 + "<<b<<"x + "<<c<<" = 0"<<endl;
	cout<<endl;
	if (a!=0)
	{
		base=((pow (b, 2))-(4*a*c));
		if (base>=0)
		{
		   raiz=pow (base, 0.5);
		   rta1=(-b+(raiz))/(2*a);
     	           rta2=(-b-(raiz))/(2*a);
    	           cout<<"La primer raiz es: "<<rta1<<endl;
    	           cout<<endl;
    	           cout<<"La segunda raiz es: "<<rta2<<endl;
		   cout<<endl;
		}
		else
		{
			img=base*(-1);
			raiz=(pow (img, 0.5))/(2*a);
			raiz2=-(pow (img, 0.5))/(2*a);
			rta1=-b/(2*a);
        	        cout<<"La primer raiz imaginaria es:  "<<rta1<<"  i "<<raiz<<endl;
			cout<<"La segunda raiz imaginaria es: "<<rta1<<"  i "<<raiz2<<endl;
        	        cout<<endl;
       	        }
	}
	else
	{
		if ((a==0) && (c==0))
		{
		  x=0/b;
		  cout<<"la raiz es: "<<x<<endl;
		}
	    	else if (b>0)
			{
	        	   lin=(-c/b);
	        	   cout<<"El valor de la raiz es de: "<<lin;
	        	   cout<<endl;
			}
		         else
				{
			            cout<<"usted solo a elegido un punto "<< c;
		        	     cout<<endl;
				}
		
	}	
}

void cubica(){

   double A,B,C,D;

   double a1,b1,c1;

   double p,q,dis;

   double x,y,z,pi;

   double ab,ac;

   double u,v,t;


   cout<<" Ecuacion cubica de la forma AX^3 + BX^2 + CX + D"<<endl;

   cout<<" Introduzca "<<endl;

   cout<<" A :"; cin>> A;

   cout<<" B :"; cin>> B;

   cout<<" C :"; cin>> C;

   cout<<" D :"; cin>> D;


   cout.precision(3); // muestra la cantidad de decimales deseada

   t = 0.0000000009; // tolerancia o margen de error

   pi = 3.141592654; // el numero pi

   a1 = B/A;

   b1 = C/A;

   c1 = D/A;

   p = b1 - (a1*a1)/3; // variable p del metodo

   q = (2*a1*a1*a1)/27 - (a1*b1)/3 + c1; // variable q del metodo

   dis = (q*q) + (4*p*p*p)/27; // discriminante

   cout<< " Ecuacion reducida : Z^3 + "<<p<< " Z + "<<q<<endl;

   cout<< " El discriminante es "<<dis<<endl;


   if ( dis >= t ) { //discriminante > 0

     ab = 0.5*(-q + sqrt(dis));

     double racub(ab); // raiz cubica 1, calculo de u

     if (ab> 0) {

	racub = exp(log(ab)/3);

     };

     if (ab== 0){

	racub = 0;

     };

     if (ab<0) {

	racub =-exp(log(-ab)/3);

     };

     ac = -0.5*(q + sqrt(dis));

	 double racub2(ac);

    if (ac> 0) {

	racub2 = exp(log(ac)/3); //raiz cubica 2, calculo de v

     };

     if (ac== 0){

	racub2 = 0;

     };

     if (ac<0) {

	racub2 =-exp(log(-ac)/3);

     };

     u = racub; // variable u del metodo de cardano

     v = racub2; // variable v del metodo de cardano

     x = u + v -(a1/3);

     y = (-0.5)*(u+v) -(a1/3);

     z = (0.5)*sqrt(3)*(u-v);

     cout<<" La ecuacion tiene una raiz real y 2 complejas "<<endl;

     cout<<" X1 :"<<x<<endl;

     cout<<" X2 :"<<y<<" + "<<z<<"i"<<endl;

     cout<<" X3 :"<<y<<" - "<<z<<"i"<<endl;

    };


   if ((dis < t and dis > -t) and (p>=-t and p<=t)){ //discriminante= p = 0

     cout<<" La ecuacion tiene una sola raiz"<<endl;

     cout<<" X :"<<-(a1/3);

	};


   if ((dis < t and dis > -t) and (p>=t or p<=-t)){ // discriminante = 0

     ab = -0.5*q;

     double racub(ab);

     if (ab> 0) {

	racub = exp(log(ab)/3); // raiz cubica de u, solo se calcula u

     };

     if (ab== 0){

	racub = 0;

     };

     if (ab<0) {

	racub =-exp(log(-ab)/3);

     };

     u = racub;

     x = 2*u -(a1/3);

     y = -u -(a1/3);

     cout<<" La ecuacion tiene raices multiples"<<endl;

     cout<<" X1    :"<<x<<endl;

     cout<<" X2=X3 :"<<y<<endl;

	};



   if ( dis <= -t ) { // discriminate < 0

    x = (2*sqrt(-p/3))*cos(acos((-q/2)*sqrt(-27/(p*p*p)))/3) -(a1/3);

    y = (2*sqrt(-p/3))*cos(acos((-q/2)*sqrt((-27)/(p*p*p)))/3 +2*pi/3) -(a1/3);

    z = (2*sqrt(-p/3))*cos(acos((-q/2)*sqrt((-27)/(p*p*p)))/3 +4*pi/3) -(a1/3);

     cout<<" La ecuacion tiene 3 raices distintas "<<endl;

     cout<<" X1 :"<<x<<endl;

     cout<<" X2 :"<<y<<endl;

     cout<<" X3 :"<<z<<endl;

    };



    getch(); // Hace una pausa al programa

} 

void eliminaciongaussiana(){

	int n, n1, i, j, k, l;

	float t, ma[10][10], x[10];



	printf("Numero de incognitas: ");

	scanf("%d", &n);

	printf("\nNumero de ecuaciones: ");

	scanf("%d", &n1);

	if(n!=n1){

		printf("No se puede ejecutar el programa. La matriz debe ser cuadratica.\n");

	}

	else{

		printf("\n");



		for(i=1;i<=n;i++){//RENGLON

			for(j=1;j<=n+1;j++){//COLUMNA

				printf("Elemento [%d][%d] de la matriz: ",i,j);

				scanf("%f", &ma[i][j]);

			}

		}



		printf("\n\nMatriz ordenada:\n\n");

		for(i=1;i<=n;i++){//IMPRIMIR MATRIZ

			for(j=1;j<=n+1;j++){

				printf("%f ", ma[i][j]);

			}

			printf("\n");

		}



		for(k=1;k<=n-1;k++){

		    int mayor=0;

                    int fmayor=k;

		    for(l=k;l<=n;l++){

                        if(fabs(ma[l][k])>mayor){

                            mayor = fabs(ma[l][k]);

                            fmayor = l;

                        }

                    }



			if(mayor == 0){

				printf("El sistema no tiene solucion.\n");

			}

			else{

            	        if(fmayor != 0){//INTERCAMBIO

                	     for(j=1;j<=n+1;j++){

                    	         t = ma[k][j];

                    	        ma[k][j] = ma[fmayor][j];

                    	        ma[fmayor][j] = t;

                	     }

            	        }

			}



                       for(i=k+1;i<=n;i++){//HACER CEROS COLUMNAS

                            t = ma[i][k] / ma[k][k];

                           for (j=k;j<=n+1;j++){

                                    ma[i][j] = ma[i][j] - t * ma[k][j];

                            }

		       }

                 }


		printf("\n\nMatriz final:\n\n");

		for(i=1;i<=n;i++){//IMPRIMIR MATRIZ FINAL

			for(j=1;j<=n+1;j++){

				printf("%f ", ma[i][j]);

			}

			printf("\n");

		}


		//ENCONTRAR RAICES

		x[n]=ma[n][n+1]/ma[n][n];



		for(i=n-1;i>=1;i--){

			float sum=0;

        	for(j=n+1;j>i;j--){

            	sum += ma[i][j] * x[j];

        	}

        	x[i] = (ma[i][n+1] - sum) / ma[i][i];

    	}

    	printf("\n\nRaices:\n\n");

		for(i=n;i>=1;i--){//IMPRIMIR RAICES

			printf("X[%d]: %f\n", i, x[i]);

		}

		printf("\n");
	}
}

void puntofijo()

{
	
printf("Programa para solucionar la funcion 2(x*x)+4x-3, con el metodo de Iteracion");

printf(" simple de Punto Fijo\n\n");

float xr1,xr2,ea,aux,fxr2;

xr1=0;

destino:

xr2=(3-(2*pow(xr1,2)))/4;

ea=fabs((xr2-xr1)/xr2);

printf("La raiz propuesta es: %f y tiene un error de: %f\n",xr2,ea);

if(ea<0.000001)

{

printf("\nLa raiz mas cercana a la solucion es: %f\n",xr2);

fxr2=(2*pow(xr2,2))+(4*xr2)-3;

printf("Y la funcion evaluada en esa raiz es: %f\n",fxr2);

}

if(ea>0.000001)

{

aux=xr2;

xr2=xr1;

xr1=aux;

goto destino;

}

}





