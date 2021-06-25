#include<vector>
#include<iostream>
#include<chrono>
#include<stdlib.h>
using namespace std;

const int count = 1000000;

class matrix_conversion{
	vector<vector<double>> transform;

	
public:
	matrix_conversion(){

	}

	vector<double> direct_conversion(vector<double>& to_be_conversed){
		using clock = std::chrono::system_clock;
using sec = std::chrono::duration<double>;
		vector<double> result(3,0);
		const auto before = clock::now();
		for(int j = 0; j<count;j++){
			
			double temp = to_be_conversed[0]/2 - (to_be_conversed[1] + to_be_conversed[2])/4;
			result[0] += temp;
			temp = to_be_conversed[0]/2 + (to_be_conversed[1] + to_be_conversed[2])/4;
			result[1] += temp;
			temp =  (to_be_conversed[2] - to_be_conversed[1])/2;
			result[2] += temp;
			
		}
		const sec duration = clock::now() - before;
		cout << "It took " << duration.count()/count << "s for direct multiplication" << endl;
		return result;
	}

	vector<double> lift_conversion(vector<double>& to_be_conversed){
		using clock = std::chrono::system_clock;
using sec = std::chrono::duration<double>;
		vector<double> result(3,0);
		const auto before = clock::now();

		for(int i = 0; i<count;i++){
			double temp1 = to_be_conversed[2] - to_be_conversed[1];

			double temp2 = to_be_conversed[1]+ (temp1/2);
			double temp3 = to_be_conversed[0] - temp2;
			double temp4 = temp2 + (temp3/2);
			
			result[0]+=(temp3/2);
			result[1]+=temp4;
			result[2]+=(temp1/2);
		}
		const sec duration = clock::now() - before;
		cout << "It took " << duration.count()/count << "s for lift_based calculation" << endl;
		return result;
	}

	void run(vector<double>& to_be_conversed){
		vector<double> r1 = direct_conversion(to_be_conversed);
		vector<double> r2 = lift_conversion(to_be_conversed);
		cout<<"result of direct_conversion():\n";
		cout<<"Cg = "<<r1[0]/count<<"  Y = "<<r1[1]/count<<" Co = "<<r1[2]/count<<endl;
		cout<<"result of lift_conversion():\n";
		cout<<"Cg = "<<r2[0]/count<<"  Y = "<<r2[1]/count<<" Co = "<<r2[2]/count<<endl;
		if(r1[0]!=r2[0]){
			cout<<"Different result for channel Cg\n";
		}
		if(r1[1]!=r2[1]){
			cout<<"Different result for channel Y\n";
		}
		if(r1[2]!=r2[2]){
			cout<<"Different result for channel Co\n";
		}



	}
};

// int main(){
// 	matrix_conversion m;
	
// 	vector<double> v = {255,55,5};
// 	m.run(v);
// }

int main(int argc, char* argv[]){
	if(argc!=4){
		cout<<"You have entered too few or too many arguments. Terminating...\n";
		return 1;
	}
	matrix_conversion m;
	
	vector<double> v(3);
	v[0]=atof((argv[1]));
	if(v[0]>255||v[0]<0){
		cout<<"The value should be in the range of [0,255]\n";
		return 1;
	}
	v[1]=atof((argv[2]));
	if(v[1]>255||v[1]<0){
		cout<<"The value should be in the range of [0,255]\n";
		return 1;
	}
	v[2]=atof((argv[3]));
	if(v[2]>255||v[2]<0){
		cout<<"The value should be in the range of [0,255]\n";
		return 1;
	}


	m.run(v);
}