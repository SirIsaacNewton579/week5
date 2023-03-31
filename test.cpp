#include<iostream>
using namespace std;
int printarr(int *p[],int m,int n){
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            cout << p[i][j];
            //cout << *(*(p+i)+j) << "\t";
        }
        cout << endl;
    }
    return 0;
}
int main(){
    int arr[3][3] = {{1,2,3},{3,4,5},{5,6,7}};
    int *p[] = {arr[0],arr[1],arr[2]};
    printarr(p,3,3);
    return 0;
}