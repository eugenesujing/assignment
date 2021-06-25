#include<vector>
#include<iostream>
#include<chrono>
#include<stdlib.h>
using namespace std;


const int count = 1;

class matrix_conversion{
	

	
public:
	matrix_conversion(){

	}

	vector<int> direct_conversion(vector<int>& to_be_conversed){
		using clock = std::chrono::system_clock;
using sec = std::chrono::duration<double>;
		vector<int> result(3,0);
		const auto before = clock::now();
		for(int j = 0; j<count;j++){
			
			int temp = (to_be_conversed[0]>>1) - ((to_be_conversed[1] + to_be_conversed[2])>>2);
			result[0] += temp;
			temp = (to_be_conversed[0]>>1) + ((to_be_conversed[1] + to_be_conversed[2])>>2);
			result[1] += temp;
			temp =  ((to_be_conversed[2] - to_be_conversed[1])>>1);
			result[2] += temp;
			
		}
		const sec duration = clock::now() - before;
		cout << "It took " << duration.count()/count << "s for direct multiplication" << endl;
		return result;
	}

	vector<int> lift_conversion(vector<int>& to_be_conversed){
		using clock = std::chrono::system_clock;
using sec = std::chrono::duration<double>;
		vector<int> result(3,0);
		const auto before = clock::now();

		for(int i = 0; i<count;i++){
			int temp1 = to_be_conversed[2] - to_be_conversed[1];

			int temp2 = to_be_conversed[1]+ (temp1>>1);
			int temp3 = to_be_conversed[0] - temp2;
			int temp4 = temp2 + (temp3>>1);
			
			result[0]+=(temp3>>1);
			result[1]+=temp4;
			result[2]+=(temp1>>1);
		}
		const sec duration = clock::now() - before;
		cout << "It took " << duration.count()/count << "s for lift_based calculation" << endl;
		return result;
	}

	void run(vector<int>& to_be_conversed){
		vector<int> r1 = direct_conversion(to_be_conversed);
		vector<int> r2 = lift_conversion(to_be_conversed);
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
	
// 	vector<int> v = {255,55,5};
// 	m.run(v);
// }

int main(int argc, char* argv[]){
	if(argc!=4){
		cout<<"You have entered too few or too many arguments. Terminating...\n";
		return 1;
	}
	matrix_conversion m;
	
	vector<int> v(3);
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