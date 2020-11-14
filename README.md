g++ -march=core2 -pthread -std=c++11 -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o SimTab-Baseline SFMT.c main.cpp

./SimTab-Baseline -d /home/liuyu/SimRace-ly-20191226/PRSim-Code-master-20191226/PRSim-Code-master/dataset/dblp-author_sorted.txt -f dblp-author -algo topk -qn 1 -k 50
