#include "MyHeader.h"
//#include "Graph.h"
//#include "NumberTheory.h"
//#include "UnionFind.h"
//#include "SegmentTree.h"
//#include "Algorithm.h"
//#include "Bit.h"

#include <vector>
#include <string>
#include <sstream>
using namespace std;

vector<string> split(const string& s, char delim) {
	vector<string> elems;
	stringstream ss(s);
	string item;
	while (getline(ss, item, delim)) {
		if (!item.empty()) {
			elems.push_back(item);
		}
	}
	return elems;
}
// ============================ Header  =================================

 
 struct TreeNode {
     int val;
     TreeNode *left;
     TreeNode *right;
     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 };
 
 class Solution {
 public:
	 int sum(TreeNode* root) {
		 int s=0;
		 if (root->left) {
			 s+=sum(root->left);
		 }
		 if (root->right) {
			 s+= sum(root->right);
		 }
		 return s;

	 }
	 TreeNode* sufficientSubset(TreeNode* root, int limit) {
		 if (root->left) {
			 if (sum(root->left) < limit)
				 root->left = NULL;
			 else
				 sufficientSubset(root->left, limit- (root->left->val));
		 }
		 if (root->right) {
			 if (sum(root->right) < limit)
				 root->right = NULL;
			 else
				 sufficientSubset(root->right, limit - (root->right->val));
		 }
		 return root;

	 }
 };


int main() {
	cin.tie(0);
	ios::sync_with_stdio(false);
	cout << fixed << setprecision(12);
	
	Solution sol;
	cout << sol.sufficientSubset("AAB"s);
	
	return 0;

}