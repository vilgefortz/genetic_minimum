#include<iostream>
#include<math.h>

using namespace std;

int main () {
	float min = 8;
	float max = 20;
	float step = (max - min)/100;
	float x = min;
	for (int i = 0; i < 100; i ++) {
		//x 8 - abs 3 - abs x 2 / +
		x+=step;
		cout << "f(" << x << ") = " << fabsf(fabsf(x-8) - 3) + (x/2) << endl;
	}
	return 0;
}

float f2 (float x) {
	return fabsf(fabsf(x-8) - 3) + (x/2)
}
