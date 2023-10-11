# The code of Estimator Accuracy Evaluation
* \accuracy\res.cpp: the row-wise sparisty estimator which is proposed by us. 
* \accuracy\meta.cpp: a statistically unbiased matrix sparsity estimation algorithm, MetaAC, which is very efficient.
* \accuracy\mnc.cpp: the state-of-the-art matrix sparsity estimation algorithm, MNC.

# The code of Efficiency Evaluation
* \efficient\resmcm.cpp: our proposed parallel SMCM algorithm.
* \efficient\rescsr.cpp: same as RESMCM except for using CSR data structure.
* \efficient\meta.cpp: same as RESMCM except for using MetaAC estimator.
* \efficient\mnc.cpp: same as RESMCM except for using MNC estimator.

# Compiling and Running
## Compiling the program
```
g++ -std=c++17 -O3 -g -fopenmp ${FILE_NAME}.cpp -o ${FILE_NAME}
```

For example:
```
g++ -std=c++17 -O3 -g -fopenmp resmcm.cpp -o resmcm
```


## Running the program:
```
nohup ${FILE_NAME} ${THREAD} ${ITERATION} ${GRAPH_FILE} ${VERTEX_FILE} ${EDGE_FILE} ${META_PATH_FILE} ${OUTPUT_FILE} 
```

For example:
```
nohup \
./resmcm \
64 \
1 \
./data/HINTrussICDE2020/DBLP/graph.txt \
./data/HINTrussICDE2020/DBLP/vertex.txt \
./data/HINTrussICDE2020/DBLP/edge.txt \
./materials/input/DBLP.txt \
./materials/output/DBLP/resmcm_ef.txt
```



## Input format
* THREAD:
The number of threads.

* ITERATION:
THE number of iteration.

* VERTEX_FILE:
Each line represents a vertex in HIN. Each line starts with the vertex_id, following by vertex type.


* EDGE_FILE:
Each line represents an edge in HIN. Each line starts with the edge_id, following by edge type.

* GRAPH_FILE:
Each line represents an adjacent array. Each line starts with the vertex_id, following by a list of neighbor_vertex_id and edge_id.

* resmcm_PATH_FILE:
Each line represents a resmcm-path.

* OUTPUT_FILE:
Each line represents a result of a resmcm-path.

    * Efficiency Evaluation: each line is consisted by "resmcm-path, nnz of output matrix, time cost"

    * Estimator Accuracy Evaluation: each line is consisted by "resmcm-path, estimated nnz"

# Other
* The five dataset used in paper are availabe from:

        DBLP: http://dblp.uni-trier.de/xml/
        DBpedia: https://wiki.dbpedia.org/Datasets
        FourSquare: https://sites.google.com/site/yangdingqi/home/foursquare-dataset
        FreeBase: http://freebase-easy.cs.uni-freiburg.de/dump/
        IMDB: https://www.imdb.com/interfaces/

* We give some 4-length meta-pathes in ```./materials/input/ ```, and the outputs of estimator accuracy evaluation and efficient evaluation of DBLP in ```./materials/out/accuracy/ ``` and ```./materials/out/efficient/ ``` respectively.
