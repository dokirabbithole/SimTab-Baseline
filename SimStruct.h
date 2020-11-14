#ifndef SIMSTRUCT_H
#define SIMSTRUCT_H

#define INT_MAX 32767

#include <vector>
#include <algorithm>
#include <queue>
#include <functional>
#include <iostream>
#include <thread>
#include <string>
#include <sstream>
#include <fstream>
#include "Graph.h"
#include "Random.h"
#include "alias.h"
#include <unordered_map>
#include <unordered_set>
#include "Timer.h"
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

/*
Make a directory
s: path of the directory
*/
int mkpath(string s, mode_t mode = 0755) {
    size_t pre = 0, pos;
    string dir;
    int mdret;
    if (s[s.size() - 1] != '/') {
        s += '/';
    }
    while ((pos = s.find_first_of('/', pre)) != string::npos) {
        dir = s.substr(0, pos++);
        pre = pos;
        if (dir.size() == 0) continue;
        if ((mdret = ::mkdir(dir.c_str(), mode)) && errno != EEXIST) {
            return mdret;
        }
    }
    return mdret;
}

// sort int in increasing order
bool comp1(const int &a, const int &b) {
	return a < b;
}

// for min heap inplemented by priority_queue
struct comp2 {
    bool operator() (int a, int b)
    {
        return a > b; // min heap
    }
};

// sort pair<int, double> in decreasing order
bool maxScoreCmp(const pair<int, double>& a, const pair<int, double>& b){
    return a.second > b.second;
}

// sort pair<int, double> in increasing order
bool minScoreCmp(const pair<int, double>& a, const pair<int, double>& b) {
    return a.second < b.second;
}

class SimStruct{
  public:
    Graph g; // Class Graph
    Random R; // Class Random
    int vert; // the number of vertices
    string filelabel; 
    double sqrtC;                                                      
    double C_value;
    double epsilon;
    double rmax;
    double nr;
    double back_nr;
    double forward_nr;
    double forward_rmax; // forward_push part's parameter 
    double backward_rmax; // backward_search's parameter
    double avg_time;
    double *H[2];//hash map,the same as U,C,UC
    int* U[2];
    int* C[1];
    int* UC[1];
    int* Parent; // 
    double *reserve;
    double *residue;
    double *newReserve;
    double *newResidue;
    bool* isInArray;

	/* for prefiltering */
	vector<pair<int, double> > estMean;
	vector<pair<int, double> > prefCB;
	double *answer;
	double *ans2;
	int *ansIdx;
	int curIdx;

	double eps_p;
	int totalSample; 
	// time
	double t_pref, t_single;

    unordered_map<int, double> dValue;

    /* data structures for MAB */
    //            query/cand   ->    level              -> w  -> \pi_l(s,w), r_l(s,w)
    unordered_map<int, unordered_map<int, unordered_map<int, double> > > adaptReserveMap, adaptResidueMap;
    vector<int> candList; // ansIdx + curIdx
    // unordered_map<int, double> mean; // answer
    // unordered_map<int, double> sum_Xi2; // ans2
    unordered_map<int, long> cntMABSample;
    long totalMABSample;
    unordered_map<int, int> budgetMap;
    unordered_map<int, double> rsumMap;
    // unordered_map<int, Alias> aliasMap;
    vector<Alias> aliasMap;
    unordered_map<int, int> nodeToAliasPos;
    double adaptPushTime, adaptDeterQueryTime, adaptRandQueryTime;

    SimStruct(string name, string file_lable, double eps, double c) {
        filelabel = file_lable;
        g.inputGraph(name);
        R = Random();
        vert = g.n;
        C_value = c;
        sqrtC = sqrt(C_value);
        epsilon = eps;
		rmax = (1-sqrtC)*(1-sqrtC)/25*epsilon;
		back_nr = (int)(20*log(vert)/epsilon/epsilon);
		nr = (int)(0.1*back_nr);

        forward_rmax = rmax * 2;
        backward_rmax = rmax;
        avg_time = 0;
        H[0] = new double[vert];
        H[1] = new double[vert];
        U[0] = new int[vert];
        U[1] = new int[vert];
        C[0] = new int[vert];
        UC[0] = new int[vert];
        Parent = new int[vert]; //
        reserve = new double[vert];
        residue = new double[vert];
        newReserve = new double[vert];
        newResidue = new double[vert];
        isInArray = new bool[vert];
        for(int i = 0; i < vert; i++){
            isInArray[i] = false;
            H[0][i] = 0;
            H[1][i] = 0;
            U[0][i] = 0;
            U[1][i] = 0;
            C[0][i] = 0;
            UC[0][i] = -1;
            Parent[i] = -1;
        }
        srand(unsigned(time(0)));
        
		// initialization
		answer = new double[vert];
		ans2 = new double[vert];
        ansIdx = new int[vert];
        for(int i = 0; i < vert; i++){
            answer[i] = 0;
			ans2[i] = 0;
            isInArray[i] = false;
        }		
		curIdx = 0;
        //tempCurIdx = 0;
        cout << "====init done!====" << endl;
    }
    ~SimStruct() {
        delete[] H[0];
        delete[] H[1];
        delete[] U[0];
        delete[] U[1];
        delete[] C[0];
        delete[] UC[0];
        delete[] Parent; //
        delete[] reserve;
        delete[] residue;
        delete[] newReserve;
        delete[] newResidue;
        delete[] isInArray;
		// for prefiltering
		delete[] answer;
		delete[] ans2;
		delete[] ansIdx;
    }

    //Estimate pi(u,w) using method FORA    
    unordered_map<int, vector<pair<int, double> > > foraMap(int u) {
        unordered_map<int, unordered_map<int, double> > answer;
        unordered_map<int, vector<pair<int, double> > > realAnswer;

        residue[u] = 1;
        vector<int> candidate_set;
        candidate_set.push_back(u);
        isInArray[u] = true;
        vector<int> new_candidate_set;
        int tempLevel = 0;

        Timer fora_timer;
        fora_timer.start();
        vector<pair<pair<int, int>, double> > aliasP;
        double rsum = 0;
        while (candidate_set.size() > 0) {
            for (int j = 0; j < candidate_set.size(); j++) {
                isInArray[candidate_set[j]] = false;
            }
            int residue_count = 0;
            for (int j = 0; j < candidate_set.size(); j++) {
                int tempNode = candidate_set[j];
                double tempR = residue[tempNode];
                newReserve[tempNode] += (1 - sqrtC) * tempR;
                int inSize = g.getInSize(tempNode);
                for (int k = 0; k < inSize; k++) {
                    int newNode = g.getInVert(tempNode, k);
                    newResidue[newNode] += residue[tempNode] * sqrtC / (double)inSize;
                    if (U[1][newNode] == 0) {
                        U[0][residue_count++] = newNode;
                        U[1][newNode] = 1;
                    }
                    if (!isInArray[newNode] && newResidue[newNode] > forward_rmax) {
                        isInArray[newNode] = true;
                        new_candidate_set.push_back(newNode);
                    }
                }
                residue[tempNode] = 0;
            }
            for (int j = 0; j < candidate_set.size(); j++) {
                int tempNode = candidate_set[j];
                if (newReserve[tempNode] > 0) {
                    answer[tempLevel][tempNode] = newReserve[tempNode];
                }
                newReserve[tempNode] = 0;
            }
            for (int j = 0; j < residue_count; j++) {
                if (newResidue[U[0][j]] <= forward_rmax) {
                    rsum += newResidue[U[0][j]];
                    aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(tempLevel + 1, U[0][j]), newResidue[U[0][j]]));
                }
                else {
                    residue[U[0][j]] = newResidue[U[0][j]];
                }
                newResidue[U[0][j]] = 0;
                U[1][U[0][j]] = 0;
            }
            candidate_set = new_candidate_set;
            new_candidate_set.clear();
            tempLevel++;
        }
        Alias alias = Alias(aliasP);
        fora_timer.end();
        cout << "forward push : " << fora_timer.timeCost << "s" << endl;
        if (rsum > 0) {
            int fora_walk_num = fora_timer.timeCost * 400000;	// 
            cout << "fora_walk_num = " << fora_walk_num << endl;
            fora_timer.start();
            double increment = rsum * (1 - sqrtC) / (double)fora_walk_num;
            // double sqrtSqrtC = sqrt(sqrtC);
            for (int j = 0; j < fora_walk_num; j++) {
                pair<int, int> tempPair = alias.generateRandom(R);
                int tempLevel = tempPair.first;
                int tempNode = tempPair.second;
                int tempCount = 0; // do not use 1 - \sqrt(1 - alpha) walk

                if (answer.find(tempLevel) != answer.end() && answer[tempLevel].find(tempNode) != answer[tempLevel].end())
                    answer[tempLevel][tempNode] += increment;
                else
                    answer[tempLevel][tempNode] = increment;
                while (R.drand() < sqrtC) { // 
                    int length = g.getInSize(tempNode);
                    if (length > 0) {
                        int r = R.generateRandom() % length;
                        tempNode = g.getInVert(tempNode, r);
                        tempCount++;
                    }
                    else {
                        break;
                    }
                    if (answer.find(tempLevel + tempCount) != answer.end() && answer[tempLevel + tempCount].find(tempNode) != answer[tempLevel + tempCount].end())
                        answer[tempLevel + tempCount][tempNode] += increment;
                    else
                        answer[tempLevel + tempCount][tempNode] = increment;
                }
            }
            fora_timer.end();
            cout << "MC : " << fora_timer.timeCost << "s" << endl;
        }
        int numFORAItems = 0, numFORAItemsPruned = 0;
        for (auto& c : answer) {
            for (auto& d : c.second) {
                realAnswer[c.first].push_back(pair<int, double>(d.first, d.second));
                numFORAItems++;
                if (d.second > forward_rmax)
                    numFORAItemsPruned++;
            }
        }
        cout << "back_nr = " << back_nr << " , numFORAItems = " << numFORAItems << " , numFORAItemsPruned = " << numFORAItemsPruned << endl;
        return realAnswer;
    }

    // estimate d(w)
    int sampleD(int nodeId, int walk_num){
        int tempCount = 0;
        for(int i = 0; i < walk_num; i++){
            int u_newNode = nodeId, v_newNode = nodeId, u_nextNode, v_nextNode;
            double meet_level = 0;
            while(R.drand() < C_value){
                int length = g.getInSize(u_newNode);
                if(length == 0)
                    break;
                int r = R.generateRandom() % length;
                u_nextNode = g.getInVert(u_newNode, r);
                length = g.getInSize(v_newNode);
                if(length == 0)
                    break;
				r = R.generateRandom() % length;
                v_nextNode = g.getInVert(v_newNode, r);
                meet_level++;
                if(u_nextNode == v_nextNode){
                    if(meet_level > 1)
                        tempCount += 1;
                    break;
                }
                u_newNode = u_nextNode;
                v_newNode = v_nextNode;
            }
        }
        return tempCount;
    }

    int sampleD2(int nodeId, int walk_num) {
        int tempCount = 0;
        for (int i = 0; i < walk_num; i++) {
            int u_newNode = nodeId, v_newNode = nodeId, u_nextNode, v_nextNode;
            while (R.drand() < C_value) {
                int length = g.getInSize(u_newNode);
                if (length == 0)
                    break;
                int r = R.generateRandom() % length;
                u_nextNode = g.getInVert(u_newNode, r);
                length = g.getInSize(v_newNode);
                if (length == 0)
                    break;
                r = R.generateRandom() % length;
                v_nextNode = g.getInVert(v_newNode, r);
                if (u_nextNode == v_nextNode) {
                    tempCount++;
                    break;
                }
                u_newNode = u_nextNode;
                v_newNode = v_nextNode;
            }
        }
        return tempCount;
    }
    
    // Calculating pi(v,w)
    void randomProbe(int w, int targetLevel, double dw, double* answer) {
        //cout << "w = " << w << " , targetLevel = " << targetLevel << endl;
        double sqrtSqrtC = sqrt(sqrtC);
        double increment = pow(sqrtSqrtC, targetLevel) * (1 - sqrtC);
        unordered_map<int, double> walkMap;

        int ind = 0;
        H[ind][w] = 1;
        int Ucount = 1;
        int Ucount1 = 0;
        int UCcount = 0;
        U[0][0] = w;
        for (int i = 0; i < targetLevel; i++) {
            //cout << "level (i) = " << i << endl;
            for (int j = 0; j < Ucount; j++) {
                int tempNode = U[ind][j];
                int outCount = g.getOutSize(tempNode);
                for (int k = 0; k < outCount; k++) {
                    int newNode = g.getOutVert(tempNode, k);
                    if (C[0][newNode] == 0) {
                        C[0][newNode] = 1;
                        UC[0][UCcount] = newNode;
                        UCcount++;
                    }
                    else {
                        C[0][newNode]++;
                    }
                }
            }

            for (int j = 0; j < UCcount; j++) {
                int tempNode = UC[0][j];
                if (R.drand() < C[0][tempNode] * sqrtSqrtC / (double)g.getInSize(tempNode)) {
                    H[1 - ind][tempNode] = 1;
                    U[1 - ind][Ucount1] = tempNode;
                    Ucount1++;
                    //cout << tempNode << endl;
                }
                C[0][UC[0][j]] = 0;
                UC[0][j] = -1;
            }
            for (int j = 0; j < Ucount; j++) {
                H[ind][U[ind][j]] = 0;
                U[ind][j] = -1;
            }
            Ucount = Ucount1;
            Ucount1 = 0;
            UCcount = 0;
            ind = 1 - ind;
            if (Ucount == 0)
                break;
        }
        
        for (int i = 0; i < Ucount; i++) {
            int tempNode = U[ind][i];
            if (walkMap.find(tempNode) == walkMap.end()) {
                walkMap[tempNode] = H[ind][tempNode] * increment;
            }
            else {
                walkMap[tempNode] += H[ind][tempNode] * increment;
            }
            U[ind][i] = 0;
            H[ind][tempNode] = 0;
        }
        Ucount = 0;

        for (auto& key : walkMap) {
            int tempNode = key.first;
            double tempSim = key.second;
            if (answer[tempNode] == 0) {
                ansIdx[curIdx] = tempNode;
                curIdx++;
            }
            double score = dw * tempSim / (1 - sqrtC) / (1 - sqrtC);
            answer[tempNode] += score;
            ans2[tempNode] += score * score;
        }
    }

    //Calculating pi(v,w)
    void randomProbeOneHopDeter(int w, int targetLevel, double dw, double* answer) {
        double sqrtSqrtC = sqrt(sqrtC);
        //double increment = pow(sqrtSqrtC, targetLevel) * (1 - sqrtC);
        unordered_map<int, double> walkMap;
        //cout << "w = " << w << " , targetLevel = " << targetLevel << endl;
        int ind = 0;
        H[ind][w] = (1 - sqrtC);
        int Ucount = 1;
        int Ucount1 = 0;
        int UCcount = 0;
        U[0][0] = w;
        for (int i = 0; i < targetLevel; i++) {
            //cout << "level (i) = " << i << endl;
            if (i < 1) {
                //cout << "deter" << endl;
                for (int j = 0; j < Ucount; j++) {
                    int tempNode = U[ind][j];
                    int outCount = g.getOutSize(tempNode);
                    for (int k = 0; k < outCount; k++) {
                        int newNode = g.getOutVert(tempNode, k);
                        //cout << "newNode = " << newNode << endl;
                        if (H[1 - ind][newNode] == 0) {
                            U[1 - ind][Ucount1] = newNode;
                            Ucount1++;
                        }
                        H[1 - ind][newNode] += sqrtC * H[ind][tempNode] / g.getInSize(newNode);
                    }
                }
                for (int j = 0; j < Ucount; j++) {
                    H[ind][U[ind][j]] = 0;
                    U[ind][j] = -1;
                }
                Ucount = Ucount1;
                Ucount1 = 0;
                UCcount = 0;
                ind = 1 - ind;
                if (Ucount == 0)
                    break;
            }
            else {
                //cout << "rand" << endl;
                for (int j = 0; j < Ucount; j++) {
                    int tempNode = U[ind][j];
                    int outCount = g.getOutSize(tempNode);
                    for (int k = 0; k < outCount; k++) {
                        int newNode = g.getOutVert(tempNode, k);
                        if (C[0][newNode] == 0) {
                            C[0][newNode] = 1;
                            Parent[newNode] = tempNode;
                            UC[0][UCcount] = newNode;
                            UCcount++;
                        }
                        else {
                            C[0][newNode]++;
                            if(R.drand() < 1.0 / C[0][newNode])
                                Parent[newNode] = tempNode;
                        }
                    }
                }

                for (int j = 0; j < UCcount; j++) {
                    int tempNode = UC[0][j];
                    if (R.drand() < C[0][tempNode] * sqrtSqrtC / (double)g.getInSize(tempNode)) {
                        //test
                        if (Parent[tempNode] == -1) {
                            cout << "error: tempNode = " << tempNode << " , Parent[tempNode] = " << Parent[tempNode] << endl;
                            exit(0);
                        } // 
                        //cout << "tempNode = " << tempNode << " , Parent[tempNode] = " << Parent[tempNode] << endl;
                        H[1 - ind][tempNode] = H[ind][Parent[tempNode]] * sqrtSqrtC;
                        U[1 - ind][Ucount1] = tempNode;
                        Ucount1++;
                    }
                    Parent[tempNode] = -1;
                    C[0][UC[0][j]] = 0;
                    UC[0][j] = -1;
                }
                for (int j = 0; j < Ucount; j++) {
                    H[ind][U[ind][j]] = 0;
                    U[ind][j] = -1;
                }
                Ucount = Ucount1;
                Ucount1 = 0;
                UCcount = 0;
                ind = 1 - ind;
                if (Ucount == 0)
                    break;
            } // end else           
        }

        for (int i = 0; i < Ucount; i++) {
            int tempNode = U[ind][i];
            if (walkMap.find(tempNode) == walkMap.end()) {
                walkMap[tempNode] = H[ind][tempNode];
            }
            else {
                walkMap[tempNode] += H[ind][tempNode];
            }
            U[ind][i] = 0;
            H[ind][tempNode] = 0;
        }
        Ucount = 0;

        for (auto& key : walkMap) {
            int tempNode = key.first;
            double tempSim = key.second;
            if (answer[tempNode] == 0) {
                ansIdx[curIdx] = tempNode;
                curIdx++;
            }
            double score = dw * tempSim / (1 - sqrtC) / (1 - sqrtC);
            answer[tempNode] += score;
            ans2[tempNode] += score * score;
        }
    }
   
    //Calculating pi(v,w) using sequential probe method
    void getRandomBackList(int w, int targetLevel, int backWalkNum, double forwardSim, double dw, double* answer) {
        unordered_map<int, double> walkMap;        
        double increment = 1 / (double)backWalkNum * pow(sqrtC, targetLevel) * (1 - sqrtC);
        for (int j = 0; j < backWalkNum; j++) {
            int ind = 0;
            H[ind][w] = 1;
            int Ucount = 1;
            int Ucount1 = 0;
            int UCcount = 0;
            U[0][0] = w;
            for (int i = 0; i < targetLevel; i++) {
                for (int j = 0; j < Ucount; j++) {
                    int tempNode = U[ind][j];
                    int outCount = g.getOutSize(tempNode);
                    double tempMaxInSize = 1 / R.drand();
                    for (int k = 0; k < outCount; k++) {
                        int newNode = g.getOutVert(tempNode, k);
                        if (g.getInSize(newNode) > tempMaxInSize) {
                            break;
                        }
                        if (H[1 - ind][newNode] == 0) {
                            U[1 - ind][Ucount1++] = newNode;
                        }
                        H[1 - ind][newNode] += H[ind][tempNode];
                    }
                }
                for (int j = 0; j < Ucount; j++) {
                    H[ind][U[ind][j]] = 0;
                    U[ind][j] = 0;
                }
                Ucount = Ucount1;
                Ucount1 = 0;
                ind = 1 - ind;
                if (Ucount == 0)
                    break;
            }
            for (int i = 0; i < Ucount; i++) {
                int tempNode = U[ind][i];
                if (walkMap.find(tempNode) == walkMap.end()) {
                    walkMap[tempNode] = H[ind][tempNode] * increment;
                }
                else {
                    walkMap[tempNode] += H[ind][tempNode] * increment;
                }
                U[ind][i] = 0;
                H[ind][tempNode] = 0;
            }
            Ucount = 0;

        }
        for (auto& key : walkMap) {
            int tempNode = key.first;
            double tempSim = key.second;
            if (answer[tempNode] == 0) {
                ansIdx[curIdx] = tempNode;
                curIdx++;
            }
            double score = dw * tempSim / (1 - sqrtC) / (1 - sqrtC);
            answer[tempNode] += score;
            ans2[tempNode] += score * score;
        }
    }
       
    /* sample-all-arms implementation */
    void singleSourceQuery(int u){
        cout << "singleSourceQuery : u = " << u << endl;
		double single_walk_num = 10 * nr;
		for(int i = 0; i < curIdx; i++){
			int nodeId = ansIdx[i];
			answer[nodeId] = 0;
			ans2[nodeId] = 0;
		}
		curIdx = 0;

        double forward_time = 0, alias_time=0, backward_time = 0, calD_time = 0, calAns_time = 0;
        unordered_map<int, vector<pair<int, double> > > forward_map = foraMap(u); // depend on forward_rmax
        vector<pair<pair<int, int>, double> > aliasP;
        double rsum = 0;
        for(auto& c: forward_map){
            int tempLevel = c.first;
            for(auto& d: c.second){
                int w = d.first;
                double forwardSim = d.second;
                rsum += forwardSim;
                aliasP.push_back(pair<pair<int, int>, double>(pair<int, int>(tempLevel,w), forwardSim));
            }
        }
        Alias alias = Alias(aliasP);
        // estimate d(w)
		unordered_map<int, int> hitMap;
        unordered_map<int, int> totalMap;
		Timer sample_d_timer;
		sample_d_timer.start();
		for(int i = 0; i < (int)single_walk_num; i++){
            pair<int, int> tempPair = alias.generateRandom(R);
			int tempLevel = tempPair.first;
            int tempNode = tempPair.second;
            int hitCount = sampleD(tempNode,1); 
			if(hitMap.find(tempNode) == hitMap.end())
                hitMap[tempNode] = 0;
            if(totalMap.find(tempNode) == totalMap.end())
                totalMap[tempNode] = 0;
            hitMap[tempNode] += hitCount;
			totalMap[tempNode] +=1;
        }
		sample_d_timer.end();
		t_single = sample_d_timer.timeCost / (int)single_walk_num;
		cout << "t_single = " << t_single << endl;

        // unordered_map<int, double> dValue;
        dValue.clear();
        int cnt_w = 0;
		for(auto& w : hitMap){
            cnt_w += 1;
			if(g.getInSize(w.first) == 0)
                dValue[w.first] = 1;
            else 
                dValue[w.first] = 1 - C_value / (double) g.getInSize(w.first) - w.second / (double)totalMap[w.first];
        }
        cout << "cal D done!" << endl;
        
		for(int i = 0; i < back_nr; i++){
			// first sample a walk from u
			int tempNode = u;
			int tempLevel = 0;
			while(R.drand() < sqrtC){
				int length = g.getInSize(tempNode);
                if(length == 0)
                    break;
                int r = R.generateRandom() % length;
                tempNode = g.getInVert(tempNode, r);
				tempLevel ++;
			}
//			pair<int, int> tempPair = alias.generateRandom(R);
//			int tempLevel = tempPair.first;
//          int tempNode = tempPair.second;
			double dw = 0;
            if(g.getInSize(tempNode) == 0)
                dw = 1;
            else if(dValue.find(tempNode) != dValue.end())
                dw = dValue[tempNode];
            else
                dw = 1 - C_value / (double) g.getInSize(tempNode);
			
			Timer pref_timer;
			pref_timer.start();
			//cout << "tempNode = " << tempNode << " , dw = " << dw << endl;
            //randomProbe(tempNode, tempLevel, dw, answer);
            //randomProbeOneHopDeter(tempNode, tempLevel, dw, answer);
            getRandomBackList(tempNode, tempLevel, 1, 1, dw, answer);
			pref_timer.end();
            t_pref += pref_timer.timeCost;
		}
		totalSample = back_nr;

		t_pref /= totalSample;
		cout << "t_pref = " << t_pref << endl;

     	answer[u] = 1.0;
    }

    double computeCB(int nodeId, double *ans2, double estMean, int totalSample, int vert){
		// invariant:
		// ans2 = \sum_{non-zero answer}{answer * answer}
		// estMean = average over all (including zero) estimations
		double failProb = 0.001; 
		double std_dev = sqrt( (ans2[nodeId] - totalSample * estMean * estMean) / totalSample );
		return std_dev * sqrt(2.0 * log(3.0 / failProb) / totalSample) + 3 * C_value * log(3.0 / failProb) / totalSample;
	}

    /*
    The prefiltering phase
    u: query node
    */
	void prefilter(int u, int k){
		cout << "prefilter : qv = " << u << endl;
        Timer prefilter_timer;
        prefilter_timer.start();
		eps_p = 1.0e-5;
		epsilon = 0.5;
		estMean.clear();
		prefCB.clear();
		
        /* iteratively apply sample-all-arms, i.e., singleSourceQuery(u) */
		do{
			//cout << "epsilon = " << epsilon << endl;
			rmax = (1 - sqrtC) * (1 - sqrtC) / 25 * epsilon;
			back_nr = (int)(20 * log(vert) / epsilon / epsilon);
			nr = (int)(0.1 * back_nr); 
			forward_rmax = rmax * 2;
			backward_rmax = rmax;
			t_pref = 0;
			t_single = 0;
			totalSample = 0;
			// invoke singleQourceQuery, should update t_pref, and t_single, totalSample, 
			clock_t local_ts = clock();
            cout << "forward_rmax = " << forward_rmax << endl;
			singleSourceQuery(u);  // baseline
			clock_t local_te = clock();
			cout << "epsilon = " << epsilon << " , time : " << (local_te - local_ts) / (double)CLOCKS_PER_SEC << endl;
			// now answer and ans2 has values
            cout << "t_pref = " << t_pref << endl;
            cout << "t_single = " << t_single << endl;
			if(t_pref > 0 && t_single > 0){
                // assume the k-th ground truth is not 0
				vector<pair<int, double> > allAnsList;
				for(int i = 0; i < curIdx; i++){
                    int nodeId = ansIdx[i];
					if(nodeId != u){
						double score = answer[nodeId] / totalSample;
						allAnsList.push_back(pair<int, double>(nodeId, score));
					}
				}
                if (allAnsList.size() >= k) {
                    sort(allAnsList.begin(), allAnsList.end(), maxScoreCmp); // in descending order, O(nlogn)
                    int kthNode = allAnsList[k - 1].first;
                    double kthScore = allAnsList[k - 1].second;
                    double kthCB = computeCB(kthNode, ans2, kthScore, totalSample, vert);
                    double cb_0 = 3 * C_value * log(1.0 / 0.0001) / totalSample;
                    cout << "kthScore = " << kthScore << " , kthCB = " << kthCB << " , kthScore - kthCB = " << kthScore - kthCB << " , cb_0 = " << cb_0 << endl;
                    if (kthScore - kthCB + eps_p > cb_0) { // can prune 0
                        vector<int> candidates;
                        for (int i = 0; i < k; i++)
                            candidates.push_back(i);
                        for (int i = k; i < allAnsList.size(); i++) {
                            double ithCB = computeCB(allAnsList[i].first, ans2, allAnsList[i].second, totalSample, vert);
                            if (kthScore - kthCB + eps_p < allAnsList[i].second + ithCB) {
                                candidates.push_back(i);
                            }
                        }
                        cout << "numCandidate = " << candidates.size() << endl;
                        if (t_pref * 100 > t_single * candidates.size()) {
                            for (int i = 0; i < candidates.size(); i++) {
                                int index = candidates[i];
                                estMean.push_back(pair<int, double>(allAnsList[index].first, allAnsList[index].second));
                                double ithCB = computeCB(allAnsList[index].first, ans2, allAnsList[index].second, totalSample, vert);
                                prefCB.push_back(pair<int, double>(allAnsList[index].first, ithCB));
                            }
                            break;
                        }
                    }
                } // otherwise keep prefiltering
			}
            else {
                cout << "t_pref = " << t_pref << " , t_single = " << t_single << endl;
                exit(0);
            }
			epsilon /= 2; 
		}while(true);
        prefilter_timer.end();
        avg_time += prefilter_timer.timeCost;
        cout << "u = " << u << " , time(sec) : " << prefilter_timer.timeCost << endl;
	}

	/* output candidate set */
	void outputCandidate(int u, int k){
		stringstream ssout;
		ssout << "candidates/" << filelabel << "/k=" << k << "/" << u << ".txt";
		ofstream fout(ssout.str());
		fout.precision(15);
		for(int i = 0; i < estMean.size(); i++){
            int nodeId = estMean[i].first;
            fout << estMean[i].first << "\t" << estMean[i].second << "\t" << totalSample << "\t" << ans2[nodeId] << "\t" << estMean.size() << endl;
		}
	}

    /*
    sample-one-arm implementation
    */
    void sampleOneArm_Naive(int queryNode, int candNode, int batchSample) {
        int tempCount = 0;
        for (int i = 0; i < batchSample; i++) {
            int u_newNode = queryNode, v_newNode = candNode, u_nextNode, v_nextNode;
            while (R.drand() < C_value) {
                int length = g.getInSize(u_newNode);
                if (length == 0)
                    break;
                int r = R.generateRandom() % length;
                u_nextNode = g.getInVert(u_newNode, r);
                length = g.getInSize(v_newNode);
                if (length == 0)
                    break;
                r = R.generateRandom() % length;
                v_nextNode = g.getInVert(v_newNode, r);

                if (u_nextNode == v_nextNode) {
                    tempCount++;
                    break;
                }
                u_newNode = u_nextNode;
                v_newNode = v_nextNode;
            }
        }
        answer[candNode] += tempCount;
        ans2[candNode] += tempCount;

        cntMABSample[queryNode] += batchSample;
        cntMABSample[candNode] += batchSample;
        totalMABSample += batchSample;
    }
    
    /* The top-k MAB phase */
    vector <pair<int, double> > UGapE(int queryNode, int k) {
        double eps_min = 1e-3;
        vector<pair<int, double> > B;
        unordered_map<int, double> UB, LB;
        vector<pair<int, double> > vecUB, vecLB;
        rsumMap.clear();
        aliasMap.clear();
        nodeToAliasPos.clear();
        // clock_t ts_ugape = clock();
        double sortTime = 0, sampleTime = 0;
        adaptPushTime = 0; adaptDeterQueryTime = 0; adaptRandQueryTime = 0;
        Timer sort_timer, sample_timer;
        // mean, cntSample, and sum_Xi2 already have values
        // clock_t tsSample = clock();
        int round = 0;
        int batch = 10;
        int batchSample = 10000; // max(10000, (int)(candList.size()));
        cout << "batchSample = " << batchSample << endl;
        // test
        //for (int i = 0; i < candList.size(); i++) {
        //    int armId = candList[i];
        //    sampleOneArm_Naive(queryNode, armId, 100);
        //} //
        while (true) {
            round++;
            sort_timer.start();
            UB.clear();
            LB.clear();
            vecUB.clear();
            vecLB.clear();
            B.clear();
            //遍历arm统一用candList
            for (int i = 0; i < candList.size(); i++) {
                int armId = candList[i];
                double ugapCB = ugapeCB(candList.size(), answer[armId], ans2[armId], cntMABSample[armId], 1.0e-4, round);
                UB[armId] = answer[armId] / cntMABSample[armId] + ugapCB;
                LB[armId] = answer[armId] / cntMABSample[armId] - ugapCB;
                vecUB.push_back(pair<int, double>(armId, answer[armId] / cntMABSample[armId] + ugapCB));
            }
            sort(vecUB.begin(), vecUB.end(), maxScoreCmp); // descreasing order
            double kthUB = vecUB[k - 1].second; // k-th largest UB
            double k1thUB = vecUB[k].second;    // k+1-th largest UB
            
            for (int i = 0; i < vecUB.size(); i++) {
                int armId = vecUB[i].first;
                if (i < k)     //(i < k)
                    B.push_back(pair<int, double>(armId, k1thUB - LB[armId]));
                else
                    B.push_back(pair<int, double>(armId, kthUB - LB[armId]));
            }
            vecUB.clear();
            sort(B.begin(), B.end(), minScoreCmp); // increasing order

            double largestB;
            for (int i = 0; i < B.size(); i++) {
                int armId = B[i].first;
                if (i < k) {
                    vecLB.push_back(pair<int, double>(armId, LB[armId]));
                    if (i == k-1)
                        largestB = B[i].second;
                }
                else
                    vecUB.push_back(pair<int, double>(armId, UB[armId]));
            }
            sort(vecLB.begin(), vecLB.end(), minScoreCmp);
            sort(vecUB.begin(), vecUB.end(), maxScoreCmp);
            sort_timer.end();
            sortTime += sort_timer.timeCost;
            
            if (largestB <= eps_min) { // B[k - 1].second
                //vector<int> upageTopk;
                vector <pair<int, double> > ugapeAnswer;
                for (int i = 0; i < vecLB.size(); i++) {
                    int armId = vecLB[i].first;
                    ugapeAnswer.push_back(pair<int, double>(armId, answer[armId] / cntMABSample[armId]));
                    cout << "armId = " << armId << " , cntMABSample[armId] = " << cntMABSample[armId] << " , est = " << answer[armId] / cntMABSample[armId] << endl;
                    //if (sameInNbr.find(armId) != sameInNbr.end()) {
                    //    vector<int>& sameInNbrList = sameInNbr[armId];
                    //    for (int j = 0; j < sameInNbrList.size(); j++)
                    //        upageTopk.push_back(sameInNbrList[j]);
                    //}
                }
                sort(ugapeAnswer.begin(), ugapeAnswer.end(), maxScoreCmp);
                //for (int i = 0; i < candList.size(); i++)
                   //delete aliasMap[i];
                
                //clock_t te_ugape = clock();
                //sampleTime = (te_ugape - tsSample) / (double)CLOCKS_PER_SEC;
                //cout << "UGapE: " << (te_ugape - ts_ugape) / (double)CLOCKS_PER_SEC << endl;
                //cout << "initTime = " << initTime << endl;
                //cout << "sampleTime = " << sampleTime << endl;
                //cout << "totalSample = " << totalSample << endl;
                //cout << "sortTime = " << sortTime << endl;
                stringstream ssout;
                ssout << "ugape/" << filelabel << "/k=" << k << "/";
                mkpath(ssout.str());
                ssout << queryNode << ".txt";
                ofstream fout(ssout.str());
                fout.precision(15);
                for (int i = 0; i < ugapeAnswer.size(); i++) {
                    fout << ugapeAnswer[i].first << "\t" << ugapeAnswer[i].second << endl;
                }
                cout << "sortTime = " << sortTime << " , sampleTime = " << sampleTime << " , adaptPushTime = " << adaptPushTime  << endl;
                cout << "adaptDeterQueryTime = " << adaptDeterQueryTime << " , adaptRandQueryTime = " << adaptRandQueryTime << endl;
                return ugapeAnswer;
            }
            // sample high
            sample_timer.start();
            for (int i = 0; i < batch; i++) {
                int armId = vecLB[i].first;
                sampleOneArm_Naive(queryNode, armId, batchSample);
                //sampleOneArm(queryNode, armId, batchSample);
            }
            // sample low
            for (int i = 0; i < batch; i++) {
                int armId = vecUB[i].first;
                sampleOneArm_Naive(queryNode, armId, batchSample);
                //sampleOneArm(queryNode, armId, batchSample);
            }
            sample_timer.end();
            sampleTime += sample_timer.timeCost;
        }
    }

    /* confidence bound */
    double ugapeCB(int numOfArms, double answer, double sum_x2, int numSample, double failProb, int round) {
        double est_mean;
        if (numSample == 0)
            est_mean = 0;
        else
            est_mean = answer / numSample;
        double estStd = sqrt((sum_x2 - numSample * est_mean * est_mean) / (numSample - 1));
        double part1 = estStd * sqrt(log(numOfArms * pow(round, 3.0) / failProb) / numSample);
        double part2 = 7.0 / 6 * log(numOfArms * pow(round, 3.0) / failProb) / (numSample - 1);
        return part1 + part2;
    }

    void readCandidatesFromPrefilter(string candpath, int qnode) {
        candList.clear();
        totalMABSample = 0;
        cntMABSample.clear();
        if (curIdx != 0) {
            for (int i = 0; i < curIdx; i++) {
                int nodeId = ansIdx[i];
                answer[nodeId] = 0;
                ans2[nodeId] = 0;
            }
            curIdx = 0;
        }
        // file format:
        // nodeId estScore totalAllArmsSample ans2 candidateSize
        stringstream ss;
        ss << candpath << "/" << qnode << ".txt";
        ifstream ifcand(ss.str());
        while (ifcand.good()) {
            int cnode;
            double estScore;
            int cntAllArmsSample;
            double ans2val;
            int candsz;
            ifcand >> cnode >> estScore >> cntAllArmsSample >> ans2val >> candsz;
            if (find(candList.begin(), candList.end(), cnode) == candList.end()) {
                candList.push_back(cnode); // equivalent to ansIdx? no must update ansIdx
                
                answer[cnode] = estScore * cntAllArmsSample;
                cntMABSample[cnode] = cntAllArmsSample;
                totalMABSample += cntAllArmsSample;
                ans2[cnode] = ans2val;

                budgetMap[cnode] = 1;
                adaptResidueMap[cnode][0][cnode] = 1;
            }
        }
        ifcand.close();
        cout << "qnode = " << qnode << " , candList.size() = " << candList.size() << endl;
        budgetMap[qnode] = 1;
        adaptResidueMap[qnode][0][qnode] = 1;

        nodeToAliasPos.clear();
        for (int i = 0; i < candList.size(); i++) {
            int nodeId = candList[i];
            ansIdx[curIdx] = nodeId;
            curIdx++;
        }
    }
};



#endif
