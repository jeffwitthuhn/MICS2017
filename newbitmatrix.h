/**
	* @file newbitmatrix.h
	* 0-1 matrix class implementation.
	* Implementations of two (0,1)-Matrix-Vector product algorithms:
	* Product carried out using Gray Codes
	* Product carried out by calculating a MST on the matrix. 
	* @authors Jeff Witthuhn
*/

#ifndef BITMATRIX
#define BITMATRIX
#include <iostream>
#include "StopWatch.h" // https://github.com/KjellKod/StopWatch
#include <unordered_map>
#include "boost/dynamic_bitset/dynamic_bitset.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/random/uniform_int_distribution.hpp"
#include <cmath>
#include <ctime>
#include <chrono> 
#include <vector>
#include <cassert>
#include "matrix.h"
#include "boost/random/uniform_real_distribution.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/properties.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/config.hpp"
#include "boost/graph/prim_minimum_spanning_tree.hpp"
#include "boost/graph/depth_first_search.hpp"
#include <functional>
#include <fstream>
#include "sequitur.hpp" //https://github.com/jsdw/cpp-sequitur

#define COUNT_OPS //to count operations instead of time

StopWatch sw; 

using namespace std;
using namespace boost; 
using namespace jw;



using DigramIndex = std::unordered_map<std::pair<SymbolWrapper,SymbolWrapper>,Symbol*>;
using RuleIndex = std::unordered_map<uint, Symbol*>;
using Value = ValueSymbol<int>;


//may be useful outside the class:
using value_type = int;
using const_value_type = const int;

const std::type_info & RuleHeadType = typeid(RuleHead);
const std::type_info & RuleTailType = typeid(RuleTail);
const std::type_info & RuleSymbolType = typeid(RuleSymbol);
const std::type_info & ValueType = typeid(Value);


// TODO: Find a better solution, perhaps figure out how to use the 
static vector<vector<int>> alist; /**< adjacency list for the MST */

// TODO: put this into the multMSTPart algorithm
static int partOpCount; /**< operation count between partitions*/
 
/** Items used for BOOST graph library algorithms (prims and DFS)*/
typedef property<edge_weight_t, int> EdgeWeightProperty;
typedef property<vertex_name_t, int> VertexNameProperty;
typedef adjacency_list<vecS, vecS, undirectedS,VertexNameProperty, EdgeWeightProperty> Graph;
typedef graph_traits<Graph>::vertex_descriptor bVertex;
typedef graph_traits<Graph>::edge_descriptor bEdge;
typedef graph_traits<Graph>::adjacency_iterator Tadjit;
typedef std::pair<Tadjit, Tadjit> viterator;

/** 
	VISITOR CLASS FOR Depth first search
	DISCOVER_VERTEX is called at every vertex visited
  this is used to convert the MST into my datastructure
*/
struct detect_loops : public boost::dfs_visitor<>
{
	
  detect_loops(int rows) {
  	alist.resize(rows);
  }
  void discover_vertex(bVertex v, const Graph& g) {
     Tadjit iter,iterend;
     tie(iter, iterend)=adjacent_vertices(v,g);
     for(;iter!=iterend;++iter) {
     	alist[v].push_back(*iter);
     }
	  return;
  }
};

/// Randomization Helpers
unsigned seed=std::chrono::system_clock::now().time_since_epoch().count();
boost::random::mt19937 gen (seed);
int random01() {
    boost::random::uniform_int_distribution<> dist(0, 1);
    return dist(gen);
}
int random0ToMax(int max) {
    boost::random::uniform_int_distribution<> dist(0, max);
    return dist(gen);
}


/**
  * \breif Main class for the algorithms BitMatrix
  *
  * Use BitMatrix class to hold implementation of Gray Code and MST
  * (0,1)-matrix-vector product algorithms
  *
  * *****************CONTROL FILE****************************
  * Must be able to do following:
  *   - Perform matrix vector product with the two new algorithms and 
  *			the naieve algorithm and report the number of addition operations
  *			required. 
  */
class BitMatrix {
	private:

	/// Used for BOOST mst and depth first search algorithms
	vector<bVertex> bVertecies;
	vector<bEdge> bEdges;
	property_map<Graph, edge_weight_t>::type HammingDistance = get(edge_weight, g);
	property_map<Graph, vertex_name_t>::type vertexname = get(vertex_name, g);
	Graph g,mst;

	bool mstCalculated;/**< true if MST calculated*/
	bool mstCalculatedPartitioned;/**< true if MST partitioned calculated*/

	std::vector<boost::dynamic_bitset<> > matrix;/**< the main datastructure */
	vector<BitMatrix> partitions; /**<partitions of matrix*/
	int rows; /**< number of rows in the matrix*/
	int columns; /**< number of columns in the matrix*/

	public: 


	/**
		default constructor
	*/
	BitMatrix() {
		rows=0;
		columns=0;
	}

	/**
		itialize size constructor
		\param rows the number of rows in the matrix
		\param columns the number of columns in the matrix
	*/
	BitMatrix(int rows, int columns) {
		this->rows=rows;
		this->columns=columns;
		mstCalculated=false;
		mstCalculatedPartitioned=false;

		matrix.resize(rows+1);
		for(int i=0;i<rows;i++)
			matrix[i].resize(columns+1);

	}

	/// Getters

	/**
		returns the specified element
		\param r row of the desired element
		\param c column of the desired element
	*/
	int getElement(int r, int c) {
		return matrix[r][c];
	}

	/**
		return the hamming distance between two rows i and j
		\param i index of the first row
		\param j index of the second row
	*/
	int getHammingDistance(int i, int j) {
		return ((matrix[i]^matrix[j]).count());
	}

	/**
		returns the number of columns in the matrix
	*/
	int numColumns() {return columns;}

	/**
		returns the number of rows in the matrix
	*/
	int numRows() {return rows;}

	/// Setters

	/**
		reinitialize the matrix
		\param rows the number of rows in the matrix
		\param columns the number of columns in the matrix
	*/
	void init(int rows, int columns) {
		this->rows=rows;
		this->columns=columns;
		this->mstCalculated=false;
		this->mstCalculatedPartitioned=false;

		matrix.resize(rows);
		for(int i=0; i<rows; i++)
			matrix[i].resize(columns);
	}

	/**
		sets all elements in the matrix to zero. 
	*/
	void setAllZero() {
		for(int i=0;i<rows;i++)
			for(int j=0;j<columns;j++)
				matrix[i][j]=1;

		matrix[1][2]=0;
	}

	/**
		sets all elements in the matrix to one. 
	*/
	void setAllOne() {
		for(int i=0;i<rows;i++)
			for(int j=0;j<columns;j++)
				matrix[i][j]=1;
	}

	/**
		sets all elements in the matrix to either 1 or zero at the desired
		density. 

		\param rate the percentage of ones in the matrix 
	*/
	void setOneDensity(int rate) {
		mstCalculated=false;
		for(int i=0;i<rows;i++)
			for(int j=0;j<columns;j++) {
				int randomNum=random0ToMax(100);
				if(randomNum<rate)
					matrix[i][j]=1;
				else 
					matrix[i][j]=0;
			}
	}

	/**
		randomize the matrix 
	*/
	void randomize() {
		mstCalculated=false;
		for(int i=0;i<rows;i++)
			for(int j=0;j<columns;j++)
				matrix[i][j]=random01();
	}

	/**
		sets a specified element to an integer value

		\param row the row index of the desired element in the matrix
		\param column the column index of the desired element in the matrix 
		\param value the value to set the desired element
	*/
	void set(int row, int column, int value) {
		matrix[row][column]=value;
		return;
	}

	/// Testing Methods

	/**
		prints a representation of the matrix to std out
	*/
	void rprint(ostream& outfile) {
		for(int i=0;i<rows;i++) {
			for(int j=0;j<columns;j++)
				outfile<<matrix[i][j];
			outfile<<std::endl;
		}
	}

	/**
		prints a graph of the hammming distances to std out
	*/
	void printHammingDistanceGraph() {
		for(int i=1;i<rows;i++) {
			for(int j=0;j<i;j++) {
				/// Magic
				std::cout<<(matrix[i]^matrix[j]).count();	
			}
			std::cout<<endl;
		}
	}

	/// Implemented Operators 

	/**
		returns the specified element
		\param r the row index of the re
	*/
	int operator() (int r, int c) {
		return matrix[r][c];
	}

	/**
		Implementation of the ordinary dense muliplication algorithm.
		Sends the count of operations to the given stream. 
		\param d the operand column vector
		\param myfile the stream to send the operation count to
	*/

	Matrix mult(Matrix d, ostream& myfile) {
		assert(columns==d.getrows());
		sw.Restart();
		int operationcount=0;
		Matrix temp(rows,1);
		double sum;
		for(int i=0;i<rows;i++) {
			sum=0;
			for(int j=0;j<columns;j++) {
				if(matrix[i][j]) {
					sum+=d.getElement(j,0);
	        #ifdef COUNT_OPS
	        	operationcount+=1;
	        #endif						
	        }
			}
			temp.set(i,0,sum);
		}
		#ifdef COUNT_OPS
		myfile << operationcount << "	";
		#else
		myfile << sw.ElapsedUs() << "	";
		#endif
		return temp;
	}

	void buildORDProgram(string filename){
		ofstream outfile;
		outfile.open(filename.c_str());
		writeHeaderORD(outfile);
		generateOperationsORD(outfile);
		writeFooterORD(outfile);
		outfile.close();
	}

	void writeHeaderORD(ostream& file){
		file 
		<< "#include <iostream>\n"
		<< "#include <fstream>\n"
		<< "#include \"StopWatch.h\"\n"
		<< "#include <string>\n"
		<< "#include <cstring>\n"
		<< "#include <vector>\n"
		<< "#include <cassert>\n"
		<< "int columns ="<<columns<< ";\n"
		<< "int rows ="<<rows<<";\n"
		<< "using namespace std;\n"
		<< "int main(int argc, char* argv[]) {\n"
		<< "	StopWatch sw;\n"
		<< "	string outFileName = argv[0];\n"
		<< "	outFileName.append(\"_out\");\n"
		<< "	ofstream outfile;\n"
		<< "	outfile.open(outFileName.c_str());\n"
		<< "	vector<double> operand(columns);\n"
		<< "	vector<double> output(rows, 0);\n\n"
		<< "	ifstream operandFile;\n"
		<< "	operandFile.open(argv[1]);\n"
		<< "	for(int i=0; i<columns; i++) {\n"
		<< "		operandFile >> operand[i];\n"
		<< "	}\n"
		<< "	// operand should now contain the operand vector\n"
		<< "	sw.Restart();\n"
		<< "	for(int i = 0; i < 10000; i++) {\n\n"
		<< "	fill(output.begin(), output.end(), 0);\n"
		<< "	// begin operations \n";

	}

	void writeFooterORD(ostream& file){
		file 
		<< "	//end operations\n\n"		
		<< "	double low = 9999999999;\n"
		<< "	double high = -9999999999;\n"
		<< "	double infnorm; \n"
		<< "	for(int j=0; j<rows; j++) {\n"
		<< "		if(output[j]>high) {\n"
		<< "			high = output[j];\n"
		<< "		}\n"
		<< "		if(output[j] < low) {\n"
		<< "			low = output[j];\n"
		<< "		}\n"
		<< "	}\n"
		<< "	if(low < 0) \n"
		<< "		low = -low;\n"
		<< "	if(high < 0)\n"
		<< "		high = -high;\n"
		<< "	if(high > low)\n"
		<< "		infnorm = high;\n"
		<< "	else\n"
		<< "		infnorm = low;\n"
		<< "	for(int j=0; j<rows; j++) {\n"
		<< "		output[j] = output[j]/infnorm;\n"
		<< "	}\n"
		<< "	operand = output; //Now use the result as new operand\n "
		<< "	}\n\n"
		<< "	outfile<<sw.ElapsedUs()/(double)CLOCKS_PER_SEC<<endl<<endl;\n"
		<< "	for(int i=0; i<rows; i++) {\n"
		<< "		outfile<<output[i]<<endl;\n"
		<< "	}\n"
		<< "	outfile.close();\n"
		<< "	return 0;\n"
		<< "}\n";
	}

	void generateOperationsORD(ostream& file) {
		for(int i=0;i<rows;i++) {
			for(int j=0;j<columns;j++) {
				if(matrix[i][j]) {
					file << "	output["<<i<<"] += operand["<<j<<"];\n";
	       }
			}
		}
	}

	/**
		Implementation of the ordinary sparse CSR muliplication algorithm.
		Sends the count of operations to the given stream. 
		\param d the operand column vector
		\param myfile the stream to send the operation count to
	*/

	Matrix multCSR(Matrix d, ostream& myfile) {
		vector<vector<int> > csr; 
		csr.resize(rows);
		for(int i =0; i<rows; i++) {
			for(int j=0; j<columns; j++) {
				if(matrix[i][j]) {
					csr[i].push_back(j);
				}
			}
		}

		assert(columns==d.getrows());
		sw.Restart();
		int operationcount=0;
		Matrix temp(rows,1);
		double sum;
		for(int i=0;i<rows;i++) {
			sum=0;
			for(int j=0;j<csr[i].size();j++) {
					sum+=d.getElement(csr[i][j],0);
					#ifdef COUNT_OPS
	        	operationcount+=1;
	        #endif				
			}
			temp.set(i,0,sum);
		}
		#ifdef COUNT_OPS
		myfile << operationcount << "	";
		#else
		myfile << sw.ElapsedUs() << "	";
		#endif
		return temp;
	}

	/**
		Helper method for the MST product algorithm
		returns XOR of two rows of the matrix 
		\param i the first row to XOR
		\param j the second row to XOR
	*/
	boost::dynamic_bitset<> XOR(int i, int j) {
		return matrix[i]^matrix[j];
	}
	
	/**
		Helper method for the MST product algorithm
		fills a Graph g with hammming distances between rows
	*/
	void fillgraph() {
		bVertecies.resize(rows);
		bEdges.resize(rows*rows);

		for(int i=0;i<rows;i++) {
			bVertecies[i]=add_vertex(g);
			vertexname[bVertecies[i]]=i;
		}

		int edgecounter=0;
		for(int i=0;i<rows;i++)
			for (int j=i+1;j<rows;j++) {
				bEdges[edgecounter]=(add_edge(bVertecies[i],bVertecies[j],g)).first;
				HammingDistance[bEdges[edgecounter]]=getHammingDistance(i,j);
				edgecounter++;
		}
	}

	/**
		calculates the MST of the hamming distance graph of the matrix. 
	*/

	void calcMST() {
		fillgraph();
		//cout<<1;
		std::vector < graph_traits < Graph >::vertex_descriptor > p(num_vertices(g));
		//cout<<2;
		prim_minimum_spanning_tree(g, &p[0]);
		//cout<<1;
		for(int i=0;i<rows;i++) {
			bVertecies[i]=add_vertex(mst);
			vertexname[bVertecies[i]]=i;

		}
		for (std::size_t i = 0; i != p.size(); ++i) {
			if(p[i]!=i)
				bEdges[i]=(add_edge(bVertecies[p[i]],bVertecies[i],mst)).first;

		}
		mstCalculated=true;
	}

	void calcMST(ostream& myfile) {
		sw.Restart();
		fillgraph();
		//cout<<1;
		std::vector < graph_traits < Graph >::vertex_descriptor > p(num_vertices(g));
		//cout<<2;
		prim_minimum_spanning_tree(g, &p[0]);
		//cout<<1;
		for(int i=0;i<rows;i++) {
			bVertecies[i]=add_vertex(mst);
			vertexname[bVertecies[i]]=i;

		}
		for (std::size_t i = 0; i != p.size(); ++i) {
			if(p[i]!=i)
				bEdges[i]=(add_edge(bVertecies[p[i]],bVertecies[i],mst)).first;

		}
		mstCalculated=true;
		myfile<<sw.ElapsedUs() << "	";
	}

	/**
		Helper method for the MST product algorithm
		performs a depth first search on the calculated MST to convert
		the MST back to my datastructure
	*/
	void calcAdjacencyList() {
		  alist.clear();
		  detect_loops vis(rows);
		  depth_first_search(mst, visitor(vis));
	}

	/**
		Computes the (0,1)-Matrix-vector product using the MST algorithm 
		sends the operation count to the given stream
		http://www.micsymposium.org/mics2016/Papers/MICS_2016_paper_12.pdf

		Algorithm: 
		- Calculate the MST of the hamming distance graph between each row.
		- Choose a starting row
		- Perform a breadth first traversal on the MST where visiting each
			node consists of calculating the corresponding element in the 
			result matrix. 
	*/
	Matrix multMST(Matrix d,ostream& myfile) {
		int operationcount=0;
		if(!mstCalculated){
			#ifdef COUNT_OPS
			calcMST();
			#else
			calcMST(myfile);
			#endif
		}
		sw.Restart();
		calcAdjacencyList();

		Matrix temp(rows,1);
		double sum=0;
		int start=minpop();

		for(int i=0;i<columns;i++) {
			if(matrix[start][i])
				sum+=d.getElement(i,0);
		}
		temp.set(start,0,sum);
		boost::dynamic_bitset<> marker;
		//TODO use queue and preallocate queues
		vector<vector<int>> tempadjlist=alist;
		vector<vector<int>> queues;
		vector<int> currentRowQueue;
		
		queues.resize(rows);
		marker.resize(rows);

		currentRowQueue.push_back(start);
		do {
			int position=currentRowQueue.back();
			currentRowQueue.pop_back();
			boost::dynamic_bitset<> differenceMatrix;
			differenceMatrix.resize(columns);
			while(!tempadjlist[position].empty()) {
				int last=tempadjlist[position].back();
				if(!marker[last]) {
						sum=temp.getElement(position,0);
						differenceMatrix=matrix[last]^matrix[position];
							for(int i=0;i<columns;i++) {
								if(differenceMatrix[i]) {
										if (matrix[position][i]) {
											sum-=d.getElement(i,0);
											#ifdef COUNT_OPS
						        	operationcount+=1;
						        	#endif				
										}
										else {
											#ifdef COUNT_OPS
						        	operationcount+=1;
						        	#endif

											sum+=d.getElement(i,0);
										}
									}
							}
						temp.set(last, 0, sum);
					currentRowQueue.push_back(last);
				}
				tempadjlist[position].pop_back();
			}
			marker[position]=1;
		}
		while(!currentRowQueue.empty());
		#ifdef COUNT_OPS
		myfile<<operationcount<<"	";
		#else 
		myfile<<sw.ElapsedUs()<<"	";
		#endif
		return temp;
	}

	/**
		finds the row with the minnimum population in the matrix. 
	*/
	int minpop() {
		int minpop=matrix[0].count();
		int minindex=0;
		int temp;
		for(int i=1;i<rows;i++) {
			temp=matrix[i].count();
			if(temp<minpop) {
				minindex=i;
				minpop=temp;
			}
		}
		return minindex;	
	}

	/**
		partition the matrix into a given number of partitions for a
		reduced required number of operations. 
		\param numberOfPartitions the number of partitions to partition the
		matrix. 
	*/
	void partition(int numberOfPartitions) {
		assert(columns%numberOfPartitions==0);
		int partitionsize=columns/numberOfPartitions;
		assert(columns%partitionsize==0);
		partitions.resize(numberOfPartitions);
		for(int i=0; i<numberOfPartitions;i++)
			partitions[i].init(rows, partitionsize);
		for(int i=0; i<numberOfPartitions;i++)
			for(int j=0;j<rows;j++)
				for (int k=partitionsize*i,m=0;m<partitionsize;m++,k++)
					partitions[i].set(j,m,getElement(j,k));
	}

	/**
		calculate the MST of each partition in the partitioned matrix
		\param numberOfPartitions the number of partitions in the matrix. 
	*/
	void calcMSTofPartitions(int numberOfPartitions) {
		sw.Restart();
		partition(numberOfPartitions);
		mstCalculatedPartitioned=true;
		for(int i=0; i<numberOfPartitions;i++) {
			partitions[i].calcMST();
		}
		cout << "calculate MST of partitions time: " << sw.ElapsedUs() << endl;
	}

	void calcMSTofPartitions(int numberOfPartitions, ostream& myfile) {
		sw.Restart();
		partition(numberOfPartitions);
		mstCalculatedPartitioned=true;
		for(int i=0; i<numberOfPartitions;i++) {
			partitions[i].calcMST();
		}
		myfile<<sw.ElapsedUs()<<"	";
	}

	/**
		Computes the (0,1)-Matrix-vector product using the MST algorithm 
		sends the operation count to the given stream. 
		Extension to the previous algorithm and first partitions the matrix
		into partitions of the given size. Then performs the algorithm on 
		each partition and combines the resultant vectors.
		http://www.micsymposium.org/mics2016/Papers/MICS_2016_paper_12.pdf 

		Algorithm: 
		- Calculate the MST of the hamming distance graph between each row.
		- Choose a starting row
		- Perform a breadth first traversal on the MST where visiting each
			node consists of calculating the corresponding element in the 
			result matrix. 
	*/
	Matrix multMSTEqPart(Matrix d,ostream& myfile, int numParts) {
		assert(columns%numParts == 0 && "partitions must divide evenly");
		return multMSTPart(d, myfile,  columns/numParts);
	}

	Matrix multMSTPart(Matrix d,ostream& myfile, int partsize) {
		int partitionsize=partsize;
		int numberOfPartitions;
		numberOfPartitions=columns/partsize;
		if(columns%partsize!=0)
			numberOfPartitions++;
		//cout<<0;
		if(!mstCalculatedPartitioned) {
			#ifdef COUNT_OPS
			sw.Restart();
			calcMSTofPartitions(numberOfPartitions);
			myfile<<sw.ElapsedUs()<<"	";
			#else
			calcMSTofPartitions(numberOfPartitions,myfile);
			#endif
		}
		//cout<<1;

		sw.Restart();
		vector<Matrix> temps;
		temps.resize(numberOfPartitions);
		for(int i=0; i<numberOfPartitions;i++) {
			temps[i].init(rows,1);
		}
		partOpCount=0;
		for(int m=0; m<numberOfPartitions;m++) {
			partitions[m].calcAdjacencyList();
				double sum=0;
				int start=partitions[m].minpop();
				for(int i=0,n=partitionsize*m;i<partitionsize;n++,i++) {


					if(partitions[m].getElement(start,i)) {
						#ifdef COUNT_OPS
						partOpCount+=1;
						#endif
						sum+=d.getElement(n,0);
					}
				}
				temps[m].set(start,0,sum);
				boost::dynamic_bitset<> marker;
				vector<vector<int>> tempadjlist=alist;
				vector<vector<int>> queues;
				vector<int> currentRowQueue;
				
				queues.resize(rows);
				marker.resize(rows);

				currentRowQueue.push_back(start);
				do{
					int position=currentRowQueue.back();
					currentRowQueue.pop_back();
					boost::dynamic_bitset<> differenceMatrix;
					differenceMatrix.resize(partitionsize);
					while(!tempadjlist[position].empty()) {
						int last=tempadjlist[position].back();
						if(!marker[last]) {
								sum=temps[m].getElement(position,0);
								differenceMatrix=partitions[m].XOR(last,position);
									for(int i=0,n=partitionsize*m;i<partitionsize;n++,i++)
									if(differenceMatrix[i]) {
											if (partitions[m].getElement(position,i)) {
												sum-=d.getElement(n,0);
												#ifdef COUNT_OPS
												partOpCount+=1;
												#endif
											}
											else {
												#ifdef COUNT_OPS
												partOpCount+=1;
												#endif

												sum+=d.getElement(n,0);
											}
										}
								temps[m].set(last, 0, sum);
							currentRowQueue.push_back(last);
						}
					tempadjlist[position].pop_back();
					}
					marker[position]=1;
				}
				while(!currentRowQueue.empty());
		}

		Matrix temp(rows,1);
		int mergenum=temp.addto(temps,numberOfPartitions);
		#ifdef COUNT_OPS
		partOpCount+=mergenum;
		#endif

		#ifdef COUNT_OPS
		myfile<<partOpCount<<"	";
		#else
		myfile<<sw.ElapsedUs()<<"	";
		#endif
		return temp;
	}

	/**
		Writes a new C++ program for the particular matrix that 
		performs the operations of the MST algorithm on a given matrix.
		The input vector is given as an argument to the program
		The output vector is written to a file. 

		This is in an attempt give the algorithm a fair timing without 
		needing to do optimizations on the algorithm implementation. 

		assumes the operand is a vector of doubles
	*/

	void buildMSTProgram(string filename, int numparts){
		ofstream outfile;
		outfile.open(filename.c_str());
		writeHeaderMST(outfile);
		generateOperationsMST(outfile, numparts);
		writeFooterMST(outfile);
		outfile.close();
	}

	void writeHeaderMST(ostream& file){
		file 
		<< "#include <iostream>\n"
		<< "#include <fstream>\n"
		<< "#include \"StopWatch.h\"\n"
		<< "#include <string>\n"
		<< "#include <cstring>\n"
		<< "#include <vector>\n"
		<< "#include <cassert>\n"
		<< "int columns ="<<columns<< ";\n"
		<< "int rows ="<<rows<<";\n"
		<< "using namespace std;\n"
		<< "int main(int argc, char* argv[]) {\n"
		<< "	StopWatch sw;\n"
		<< "	string outFileName = argv[0];\n"
		<< "	outFileName.append(\"_out\");\n"
		<< "	ofstream outfile;\n"
		<< "	outfile.open(outFileName.c_str());\n"
		<< "	vector<double> operand(columns);\n"
		<< "	vector<double> temp(rows, 0);\n"
		<< "	vector<double> output(rows, 0);\n\n"
		<< "	ifstream operandFile;\n"
		<< "	operandFile.open(argv[1]);\n"
		<< "	for(int i=0; i<columns; i++) {\n"
		<< "		operandFile >> operand[i];\n"
		<< "	}\n"
		<< "	// operand should now contain the operand vector\n"
		<< "	sw.Restart();\n"
		<< "	for(int i = 0; i < 10000; i++) {\n\n"
		<< "	fill(temp.begin(), temp.end(), 0);\n"
		<< "	fill(output.begin(), output.end(), 0);\n"
		<< "	// begin operations \n";

	}

	void writeFooterMST(ostream& file){
		file 
		<< "	//end operations\n\n"		
		<< "	double low = 9999999999;\n"
		<< "	double high = -9999999999;\n"
		<< "	double infnorm; \n"
		<< "	for(int j=0; j<rows; j++) {\n"
		<< "		if(output[j]>high) {\n"
		<< "			high = output[j];\n"
		<< "		}\n"
		<< "		if(output[j] < low) {\n"
		<< "			low = output[j];\n"
		<< "		}\n"
		<< "	}\n"
		<< "	if(low < 0) \n"
		<< "		low = -low;\n"
		<< "	if(high < 0)\n"
		<< "		high = -high;\n"
		<< "	if(high > low)\n"
		<< "		infnorm = high;\n"
		<< "	else\n"
		<< "		infnorm = low;\n"
		<< "	for(int j=0; j<rows; j++) {\n"
		<< "		output[j] = output[j]/infnorm;\n"
		<< "	}\n"
		<< "	operand = output; //Now use the result as new operand\n "
		<< "	}\n\n"
		<< "	outfile<<sw.ElapsedUs()/(double)CLOCKS_PER_SEC<<endl<<endl;\n"
		<< "	for(int i=0; i<rows; i++) {\n"
		<< "		outfile<<output[i]<<endl;\n"
		<< "	}\n"
		<< "	outfile.close();\n"
		<< "	return 0;\n"
		<< "}\n";
	}
		

	void writeOperation(ostream& file, int indexIn, int indexOut, int op) {
		if(op == 0) {
			file << "	temp["<<indexOut<<"] -= operand["<<indexIn<<"];\n";
		}
		else if(op == 1){
			file << "	temp["<<indexOut<<"] += operand["<<indexIn<<"];\n";
		}
		else if(op == 2) {
			file << "	temp["<<indexOut<<"] = temp["<<indexIn<<"];\n";
		}
		else {
			file 
			<< "	for(int i=0; i<rows; i++) {\n"
			<< "			output[i] += temp[i];\n"
			<< "			temp[i] = 0;\n"
			<< "		}\n";
			;
		}
	}

	void generateOperationsMST(ostream& file, int numParts) {
		assert(columns%numParts == 0 && "partitions must divide evenly");
		int partsize = columns/numParts;
		int partitionsize=partsize;
		int numberOfPartitions;
		numberOfPartitions=columns/partsize;
		if(columns%partsize!=0)
			numberOfPartitions++;

		if(!mstCalculatedPartitioned) {
			calcMSTofPartitions(numberOfPartitions);
		}

		sw.Restart();
		vector<Matrix> temps;
		temps.resize(numberOfPartitions);
		for(int i=0; i<numberOfPartitions;i++) {
			temps[i].init(rows,1);
		}
		partOpCount=0;
		for(int m=0; m<numberOfPartitions;m++) {
			partitions[m].calcAdjacencyList();
				double sum=0;
				int start=partitions[m].minpop();
				for(int i=0,n=partitionsize*m;i<partitionsize;n++,i++) {
					if(partitions[m].getElement(start,i)) {
						writeOperation(file,n, start, 1);
					}
				}
				boost::dynamic_bitset<> marker;
				vector<vector<int>> tempadjlist=alist;
				vector<vector<int>> queues;
				vector<int> currentRowQueue;
				
				queues.resize(rows);
				marker.resize(rows);

				currentRowQueue.push_back(start);
				do{
					int position=currentRowQueue.back();
					currentRowQueue.pop_back();
					boost::dynamic_bitset<> differenceMatrix;
					differenceMatrix.resize(partitionsize);
					while(!tempadjlist[position].empty()) {
						int last=tempadjlist[position].back();
						if(!marker[last]) {
								sum=temps[m].getElement(position,0);
								writeOperation(file,position, last, 2);
								differenceMatrix=partitions[m].XOR(last,position);
									for(int i=0,n=partitionsize*m;i<partitionsize;n++,i++)
										if(differenceMatrix[i]) {
												if (partitions[m].getElement(position,i)) {
													writeOperation(file,n, last, 0);
												}
												else {
													writeOperation(file,n, last, 1);
												}
											}
							currentRowQueue.push_back(last);
						}
					tempadjlist[position].pop_back();
					}
					marker[position]=1;
				}
				while(!currentRowQueue.empty());
				writeOperation(file,0,0,3);
		}
	}




	/// Gray Code Algorithm Methods
	int grayRows; /**< number fo rows in the graycode*/
	int numParts; /**< number of partitions in the gray code*/
	std::vector<boost::dynamic_bitset<> > gray; /**< the gray code*/

	/**
		builds/initializes reflective gray code 
		\param numberOfPartitions the number of partitions in the matrix. 
	*/
	void initGrayCode() {
		grayRows = 1<<columns;
		gray.resize(grayRows, boost::dynamic_bitset<>(columns));
		for(int i=0; i<grayRows; i++) {
			for(int j=0; j<columns;j++) {
				gray[i][j]=(i>>j)&1;
			}
			gray[i]=gray[i]^(gray[i]>>1);
		}
	}

	/**
		Uses a Gray Code to comptuer the (0,1)-matrix-vector product by 
		first computing the Gray Code-vector product which can be done 
		quickly since traversing the gray code each row has a hamming 
		distance of one. Then the matrix-vector product is computed by 
		looking up the corresponding rows in the gray code product result.

		http://www.micsymposium.org/mics2016/Papers/MICS_2016_paper_12.pdf

		the matrix is partitioned into partitions of the given size and 
		then the algorithm is done on each partition and the resultant 
		vectors are combined. 
	*/
	Matrix multGrayParts(Matrix d, ostream& myfile, int partSize) {
		int start;
		int opcount=0; 
		assert(columns%partSize==0&&"columns must be a multiple of partSize"); 
		Matrix output(rows, 1);
		output.setzero();
		numParts = columns / partSize;
		int grayPartRows = 1<<partSize;
		std::vector<boost::dynamic_bitset<> > grayPart(grayPartRows, boost::dynamic_bitset<>(partSize));
		for(int i=0; i<grayPartRows; i++) {
			for(int j=0; j<partSize;j++) {
				grayPart[i][j]=(i>>j)&1;
			}
			grayPart[i]=grayPart[i]^(grayPart[i]>>1);
		}
		std::vector<Matrix> scratch(numParts, Matrix(grayPartRows,1));
		boost::dynamic_bitset<> differset;
		for(int k=0; k<numParts; k++) {
			start = k * partSize;
			for(int i=0; i<grayPartRows; i++) {

				if(i==0)
					scratch[k].set(0,0,0);
				else{
					differset=grayPart[i]^grayPart[i-1];
					for(int j=0; j<partSize;j++) {
					if(differset[j]) {
						opcount+=1;
						if(grayPart[i][j])
							scratch[k].set(i,0,scratch[k].getElement(i-1,0)+d.getElement(start+j,0));
						else 
							scratch[k].set(i,0,scratch[k].getElement(i-1,0)-d.getElement(start+j,0));
						break;
						}
					}
				}
			}
		}

		boost::dynamic_bitset<> G(partSize),index, entireRow;
		int index_int;
		for(int k=0; k<numParts;k++) {
			for(int i=0; i<rows; i++) {
				entireRow=matrix[i];
				start = k*partSize;
				for(int j=0, jk=start; j<partSize; j++, jk++)
					G[j]=entireRow[jk];
				index=G;
				for(int j=1; j<partSize+i; j++) {
					index=index^(G>>j);
				}
				index_int=0;
				for(int j=0; j<index.size(); j++) {
					index_int=index_int|(index[j]<<j);
				}
				if(scratch[k].getElement(index_int,0)) {
					opcount++;
					output.set(i,0,scratch[k].getElement(index_int,0)+output.getElement(i,0));
				}
			}
		}

		myfile<<opcount<<"	";
		return output;
	}
	/**
		Uses a Gray Code to comptuer the (0,1)-matrix-vector product by 
		first computing the Gray Code-vector product which can be done 
		quickly since traversing the gray code each row has a hamming 
		distance of one. Then the matrix-vector product is computed by 
		looking up the corresponding rows in the gray code product result.

		http://www.micsymposium.org/mics2016/Papers/MICS_2016_paper_12.pdf
	*/
	Matrix multGray(Matrix d,ostream& myfile) {
		boost::dynamic_bitset<> differset;
		Matrix output(rows, 1);
		int opcount=0;
		initGrayCode();
		Matrix scratch(grayRows,1);
		for(int i=0; i<grayRows; i++) {	
			if(i==0)
				scratch.set(0,0,0);
			else{
				differset=gray[i]^gray[i-1];
				for(int j=0; j<columns;j++) {
				if(differset[j]) {
					opcount++;
					if(gray[i][j])
						scratch.set(i,0,scratch.getElement(i-1,0)+d.getElement(j,0));
					else 
						scratch.set(i,0,scratch.getElement(i-1,0)-d.getElement(j,0));
					break;
					}
				}
			}
		}

		boost::dynamic_bitset<> G,index;
		int index_int;
		
		for(int i=0; i<rows; i++) {
			G=matrix[i];
			index=G;
			for(int j=1; j<columns; j++) {
				index=index^(G>>j);
			}
			index_int=0;
			for(int j=0; j<index.size(); j++) {
				index_int=index_int|(index[j]<<j);
			}
			output.set(i,0,scratch.getElement(index_int,0));
		}
		myfile<<opcount<<"	";
		return output;

	}

	///SEQUITER algorithm
	/**
		Matrix vector product via compression using SEQUITER algorithm.
		method described in:
		http://www.micsymposium.org/mics_2005/papers/paper60.pdf
		
		build grammar using a compressed sparse row representation,
		only keeping track os 1's
	*/
	void buildMSTandSeq() {
		cout<<"starting precalculations\n";
		calcMST();
		cout<<"MST CALCULATED, starting seq\n";
		buildSeq();
		cout<<"SEQ BUILT\n";
	}

	Sequitur<int> seq;
	std::unordered_map<int,int> ruleValueMap; 
	std::unordered_map<uint, Symbol*> rules;
	bool seqBuilt = false;
  int seqopcount;
  /**
		generate the rules for the grammar
  */
	void buildSeq() {
		sw.Restart();
		for(int i=0; i<rows; i++) {
			for(int j=0; j<columns; j++) {
				if(matrix[i][j])
					seq.push_back(j);
			}
			seq.push_back(-1-i);
		}
		rules = seq.getRules();
		seqBuilt = true;
		seqopcount=0; 
		cout<<"build seq time: " << sw.ElapsedUs() << endl;
	}

	void buildSeq(ostream& myfile) {
		sw.Restart();
		for(int i=0; i<rows; i++) {
			for(int j=0; j<columns; j++) {
				if(matrix[i][j])
					seq.push_back(j);
			}
			seq.push_back(-1-i);
		}
		rules = seq.getRules();
		seqBuilt = true;
		seqopcount=0; 
		myfile<<sw.ElapsedUs()<<"	";
	}


	/**
	traverse the grammar recursively, computing the matrix-vector product. 
	*/
		void buildSEQProgram(string filename){
		ofstream outfile;
		outfile.open(filename.c_str());
		writeHeaderSEQ(outfile);
		generateOperationsSEQ(outfile);
		writeFooterSEQ(outfile);
		outfile.close();
	}

	void writeHeaderSEQ(ostream& file){
		buildSeq();
		int size = countRules() + 1;
		file 
		<< "#include <iostream>\n"
		<< "#include <fstream>\n"
		<< "#include \"StopWatch.h\"\n"
		<< "#include <string>\n"
		<< "#include <cstring>\n"
		<< "#include <vector>\n"
		<< "#include <cassert>\n"
		<< "int columns ="<<columns<< ";\n"
		<< "int rows ="<<rows<<";\n"
		<< "int ruleSize ="<<size<<";\n"
		<< "using namespace std;\n"
		<< "int main(int argc, char* argv[]) {\n"
		<< "	StopWatch sw;\n"
		<< "	string outFileName = argv[0];\n"
		<< "	outFileName.append(\"_out\");\n"
		<< "	ofstream outfile;\n"
		<< "	outfile.open(outFileName.c_str());\n"
		<< "	vector<double> operand(columns);\n"
		<< "	vector<double> ruleValues(ruleSize, 0);\n"
		<< "	vector<double> output(rows, 0);\n"
		<< "	ifstream operandFile;\n"
		<< "	operandFile.open(argv[1]);\n"
		<< "	for(int i=0; i<columns; i++) {\n"
		<< "		operandFile >> operand[i];\n"
		<< "	}\n"
		<< "	// operand should now contain the operand vector\n"
		<< "	sw.Restart();\n"
		<< "	for(int i = 0; i < 10000; i++) {\n\n"
		<< "	fill(ruleValues.begin(), ruleValues.end(), 0);\n"
		<< "	fill(output.begin(), output.end(), 0);\n"
		<< "	// begin operations \n";

	}

	void writeFooterSEQ(ostream& file){
		file 
		<< "	//end operations\n\n"		
		<< "	double low = 9999999999;\n"
		<< "	double high = -9999999999;\n"
		<< "	double infnorm; \n"
		<< "	for(int j=0; j<rows; j++) {\n"
		<< "		if(output[j]>high) {\n"
		<< "			high = output[j];\n"
		<< "		}\n"
		<< "		if(output[j] < low) {\n"
		<< "			low = output[j];\n"
		<< "		}\n"
		<< "	}\n"
		<< "	if(low < 0) \n"
		<< "		low = -low;\n"
		<< "	if(high < 0)\n"
		<< "		high = -high;\n"
		<< "	if(high > low)\n"
		<< "		infnorm = high;\n"
		<< "	else\n"
		<< "		infnorm = low;\n"
		<< "	for(int j=0; j<rows; j++) {\n"
		<< "		output[j] = output[j]/infnorm;\n"
		<< "	}\n"
		<< "	operand = output; //Now use the result as new operand\n "
		<< "	}\n\n"
		<< "	outfile<<sw.ElapsedUs()/(double)CLOCKS_PER_SEC<<endl<<endl;\n"
		<< "	for(int i=0; i<rows; i++) {\n"
		<< "		outfile<<output[i]<<endl;\n"
		<< "	}\n"
		<< "	outfile.close();\n"
		<< "	return 0;\n"
		<< "}\n";
	}

	void writeOperationSEQ(ostream& file, int indexIn, int indexOut, int op) {
		if(op == 0) {
			file << "	output["<<indexOut<<"] += operand["<<indexIn<<"];\n";
		}
		else if(op == 1) {
			file << "	output["<<indexOut<<"] += ruleValues["<<indexIn<<"];\n";

		}
		else if(op == 2) {
			file << "	ruleValues["<<indexOut<<"] += ruleValues["<<indexIn<<"];\n";
		}
		else if(op == 3) {
			file << "	ruleValues["<<indexOut<<"] += operand["<<indexIn<<"];\n";
		}
		
		
	}

	void generateOperationsSEQ(ostream& file) {
		ruleValueMap.clear();
		unsigned int number = 4; 
		sw.Restart();
		Matrix ans(columns, 1);
		int i =0; 
		auto current = rules[0];
		while (i < rows) {
			bool done = false;
			double sum = 0;
			while(current && !done) {
		    if(typeid(*current) == seq.RuleSymbolType) {
		        auto r = static_cast<RuleSymbol*>(current);
		       	computeRuleValueProgram(file, rules[r->getID()], r->getID());
		       	writeOperationSEQ(file,r->getID(), i,1);
		      }
		    else if(typeid(*current) == seq.ValueType) {
		        auto value = static_cast<ValueSymbol<int>*>(current);
		        if(value->getValue() < 0) {
		        	done = true;
		        }
		        else {
		       		writeOperationSEQ(file,value->getValue(), i,0);
		        }
		      }
		   	current = current -> next();
		  }
		  ans.set(i,0,sum);
		  if(i%(rows/10) == 0)
		  	cout<<'|';
			i++;
		}
		cout<<endl;
	}

	void computeRuleValueProgram(ostream& file,Symbol* rule, unsigned int index) {
		if (ruleValueMap.find(index) != ruleValueMap.end()) {
			return;
		}

		double sum = 0; 
		auto current = rule; 
		
	  while(current) {
	    if(typeid(*current) == seq.RuleSymbolType) {
	        auto r = static_cast<RuleSymbol*>(current);
	        computeRuleValueProgram(file, rules[r->getID()], r->getID());
	        writeOperationSEQ(file,r->getID(), index, 2);
         
		      }
	    else if(typeid(*current) == seq.ValueType) {
	        auto value = static_cast<ValueSymbol<int>*>(current);
	        if(value->getValue() < 0) {
	        	cout<<"ERROR shouldnt find negative indeces: ";
	        	cout<<value->getValue() << endl;
	        	printRules();
	          auto r = static_cast<RuleSymbol*>(rule);

	        	printRule(r->getID());
	        }
	        else {
	         	//sum += d.getElement(value->getValue(),0);
	        	writeOperationSEQ(file,value->getValue(), index, 3);

		       }

	      }
	   	current = current -> next();
	  }
	  ruleValueMap.emplace(index, sum);
	}


	Matrix multSequiter(Matrix d,ostream& myfile) {
		ruleValueMap.clear();
		if(!seqBuilt){
			#ifdef COUNT_OPS
			buildSeq();
			#else
			buildSeq(myfile);
			#endif
		}
		unsigned int number = 4; 
		sw.Restart();
		Matrix ans(columns, 1);
		int i =0; 
		auto current = rules[0];
		while (i < rows) {
			bool done = false;
			double sum = 0;
			while(current && !done) {
		    if(typeid(*current) == seq.RuleSymbolType) {
		        auto r = static_cast<RuleSymbol*>(current);
		        sum+=computeRuleValue(d, rules[r->getID()], r->getID());
		        #ifdef COUNT_OPS
		        	seqopcount+=1;
		        #endif
		        }
		    else if(typeid(*current) == seq.ValueType) {
		        auto value = static_cast<ValueSymbol<int>*>(current);
		        if(value->getValue() < 0) {
		        	done = true;
		        }
		        else {
		         sum += d.getElement(value->getValue(),0);
		        #ifdef COUNT_OPS
		        	seqopcount+=1;
		        #endif		        
		        }
		      }
		   	current = current -> next();
		  }
		  ans.set(i,0,sum);
		  if(i%(rows/10) == 0)
		  	cout<<'|';
			i++;
		}
		cout<<endl;
		#ifdef COUNT_OPS
		myfile<<seqopcount<<"	";
		#else
		myfile<<sw.ElapsedUs()<<"	";
		#endif
		return ans;

	}
	
	int computeRuleValue(Matrix d, Symbol* rule, unsigned int index) {
		if (ruleValueMap.find(index) != ruleValueMap.end()) {
			return ruleValueMap[index];
		}

		double sum = 0; 
		auto current = rule; 
	  while(current) {
	    if(typeid(*current) == seq.RuleSymbolType) {
	        auto r = static_cast<RuleSymbol*>(current);
	        sum+=computeRuleValue(d, rules[r->getID()], r->getID());
		        #ifdef COUNT_OPS
		        	seqopcount+=1;
		        #endif	       
		      }
	    else if(typeid(*current) == seq.ValueType) {
	        auto value = static_cast<ValueSymbol<int>*>(current);
	        if(value->getValue() < 0) {
	        	cout<<"ERROR shouldnt use rule 0\n";
	        }
	        else {
	         sum += d.getElement(value->getValue(),0);
		        #ifdef COUNT_OPS
		        	seqopcount+=1;
		        #endif	       

		         }

	      }
	   	current = current -> next();
	  }
	  ruleValueMap.emplace(index, sum);
		return sum;
	}

	void printRule(unsigned int index) {
		auto rule = rules[index];
		  
		  auto current = rule;
		  while(current) {
		      if(typeid(*current) == seq.RuleHeadType) {
		          auto head = static_cast<RuleHead*>(current);
		          cout << "R" <<head->getID()<<"("<<head->getCount()<<"): ";
		          }
		      else if(typeid(*current) == seq.RuleSymbolType) {
		          auto rule = static_cast<RuleSymbol*>(current);
		          cout << "R" << rule->getID() << " ";
		          }
		      else if(typeid(*current) == seq.ValueType) {
		          auto value = static_cast<ValueSymbol<int>*>(current);
		          if(value->getValue() < 0) {
		          	cout<<"|" << " ";
		          }
		          else {
		         	 cout << "C" << value->getValue() << " ";
		          }
		        }
		      else if (typeid(*current) == seq.RuleTailType) {
		          cout << endl;
		          }
		      current = current -> next();
		      }
	}

	/**
	 return number of rules in the grammar
	*/
	int countRules() {
		int i=0;
		for(auto rule: rules) {
			i++;
		}
		return i;
	}
	/**
		helper function to print the rules in the grammar. 
	*/
	void printRules()  {
		for(auto rule : rules)
		  {
		  //rule is of type <unsigned,Symbol*>
		  
		  //print rule ID:
		  cout << "rule ID: " << rule.first << endl;
		  
		  //get first rule symbol:
		  auto current = rule.second;
		  while(current)
		      {
		      if(typeid(*current) == seq.RuleHeadType)
		          {
		          auto head = static_cast<RuleHead*>(current);
		          cout << "Rule Head" <<endl;
		          cout << "- count: " << head->getCount() << endl;
		          cout << "- ID: " << head->getID() << endl;
		          }
		      else if(typeid(*current) == seq.RuleSymbolType)
		          {
		          auto rule = static_cast<RuleSymbol*>(current);
		          cout << "Pointer to rule: " << rule->getID() << endl;
		          }
		      else if(typeid(*current) == seq.ValueType)
		          {
		          auto value = static_cast<ValueSymbol<int>*>(current);
		          cout << "Value: " << value->getValue() << endl;
		          }
		      else
		          {
		          cout << "Rule Tail" << endl;
		          }
		          current = current -> next();
		      }
		  }
	}


	///SEQUITER algorithm  partitioned
	/**
		Matrix vector product via compression using SEQUITER algorithm.
		method described in:
		http://www.micsymposium.org/mics_2005/papers/paper60.pdf
	
		partitiones the problem and comptutes sub products, 
		see MST and gray code partitioned. 
	*/

	vector<Sequitur<int>> seqParts;
	vector<std::unordered_map<int,int>> ruleValueMapParts; 
	vector<std::unordered_map<uint, Symbol*>> ruleParts;
	bool partseqBuilt = false;
	int numPartsSeq; 
  int partseqopcount;

  /**
		generate the rules for the grammar
  */
	void buildSeqParts(int partSize) {
		assert(columns%partSize == 0 && "assert columns%partSize == 0");
		numPartsSeq= columns/partSize; 
		seqParts.resize(numPartsSeq);
		ruleValueMapParts.resize(numPartsSeq);
		ruleParts.resize(numPartsSeq);
		sw.Restart();
		for(int k=0; k<numPartsSeq; k++) {
			for(int i=0; i<rows; i++) {
				for(int j=k*partSize; j<(k+1)*partSize; j++) {
					if(matrix[i][j])
						seqParts[k].push_back(j);
				}
				seqParts[k].push_back(-1-i);
			}
			ruleParts[k] = seqParts[k].getRules();
		}
		
		partseqBuilt = true;
		partseqopcount=0; 
		cout<<"build seqParts time: " << sw.ElapsedUs() << endl;
	}

	void buildSeqParts(int partSize, ostream& myfile) {
		assert(columns%partSize == 0 && "assert columns%partSize == 0");
		numPartsSeq= columns/partSize; 
		seqParts.resize(numPartsSeq);
		ruleValueMapParts.resize(numPartsSeq);
		ruleParts.resize(numPartsSeq);
		sw.Restart();
		for(int k=0; k<numPartsSeq; k++) {
			for(int i=0; i<rows; i++) {
				for(int j=k*partSize; j<(k+1)*partSize; j++) {
					if(matrix[i][j])
						seqParts[k].push_back(j);
				}
				seqParts[k].push_back(-1-i);
			}
			ruleParts[k] = seqParts[k].getRules();
		}
		
		partseqBuilt = true;
		partseqopcount=0; 
		myfile<<sw.ElapsedUs() <<"	";
	}
	/**
	traverse the grammar recursively, computing the matrix-vector product. 
	TODO: instead of traversing linked list, convert to an array first in build.
	*/
	Matrix multSequiterEqParts(Matrix d, ostream& myfile, int numParts) {
		assert(columns%numParts == 0 && "partitions must divide evenly");
		return multSequiterParts(d, myfile,  columns/numParts);
	}
	Matrix multSequiterParts(Matrix d, ostream& myfile, int partSize) {
		ruleValueMap.clear();
		if(!partseqBuilt){
			#ifdef COUNT_OPS
			sw.Restart();
			buildSeqParts(partSize);
			myfile<<sw.ElapsedUs()<<"	";
			#else
			buildSeqParts(partSize, myfile);
			#endif	
		}
		unsigned int number = 4; 
		cout<<"numparts: " << numPartsSeq << endl;
		sw.Restart();
		Matrix ans(columns, 1);
		for(int k=0; k<numPartsSeq; k++){
			#ifdef COUNT_OPS
			cout<<" [k:" << k << "] |";
			#endif
			int i =0; 
			auto current = ruleParts[k][0];
			while (i < rows) {
				bool done = false;
				double sum = 0;
				while(current && !done) {
			    if(typeid(*current) == seq.RuleSymbolType) {
			        auto r = static_cast<RuleSymbol*>(current);
			        sum+=computeRuleValueParts(d, ruleParts[k][r->getID()], r->getID(),k);
			        #ifdef COUNT_OPS
			        	partseqopcount+=1;
			        #endif			        
			        }
			    else if(typeid(*current) == seq.ValueType) {
			        auto value = static_cast<ValueSymbol<int>*>(current);
			        if(value->getValue() < 0) {
			        	done = true;
			        }
			        else {
			         sum += d.getElement(value->getValue(),0);
				        #ifdef COUNT_OPS
				        	partseqopcount+=1;
				        #endif					      
				        }
			      }
			   	current = current -> next();
			  }
			  ans.set(i,0,ans.getElement(i,0)+ sum);
			#ifdef COUNT_OPS
			  if(i%(rows/10) == 0)
			  	cout<<'|';
			 #endif
				i++;
			}
		}
		#ifdef COUNT_OPS
		cout<<endl;
		myfile<<partseqopcount<<"	";
		#else
		myfile<<sw.ElapsedUs()<<"	";
		#endif
		return ans;
	}

	int computeRuleValueParts(Matrix d, Symbol* rule, unsigned int index, int k) {
		if (ruleValueMapParts[k].find(index) != ruleValueMapParts[k].end()) {
			return ruleValueMapParts[k][index];
		}

		double sum = 0; 
		auto current = rule; 
	  while(current) {
	    if(typeid(*current) == seq.RuleSymbolType) {
	        auto r = static_cast<RuleSymbol*>(current);
	        sum+=computeRuleValueParts(d, ruleParts[k][r->getID()], r->getID(),k);
		        #ifdef COUNT_OPS
		        	partseqopcount+=1;
		        #endif			        
		      }
	    else if(typeid(*current) == seq.ValueType) {
	        auto value = static_cast<ValueSymbol<int>*>(current);
	        if(value->getValue() < 0) {
	        	cout<<"ERROR shouldnt use rule 0\n";
	        }
	        else {
	         sum += d.getElement(value->getValue(),0);
		        #ifdef COUNT_OPS
		        	partseqopcount+=1;
		        #endif		
	        }
	      }
	   	current = current -> next();
	  }
	  ruleValueMapParts[k].emplace(index, sum);
		return sum;
	}


};

#endif
