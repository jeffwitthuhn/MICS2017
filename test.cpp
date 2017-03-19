/**
  * @file test.cppl
  * testing the Matrix algorithm implementations
  * @authors Jeff Witthuhn
*/

#include "newbitmatrix.h"
#include "sequitur.hpp" //https://github.com/jsdw/cpp-sequitur
#include <string>
#include <fstream>
#include <chrono>
#include <unordered_map>
#include "StopWatch.h"
#include <stdlib.h>
#include <cmath>
#include <time.h>

//compile with: g++ test.cpp -IC:\boost\boost_1_59_0 -std=c++11
using namespace std;

/**
	DATA COLLECTION Naieve VS SEQUITUR VS MST algorthms
	Will look at the operation count when varying, 
	 num rows to num column ratio, sparsity, number of elements.

	default values (that are kept constant unless testing that):
	ROWS/COLUMNS==1
	NUMBER OF ELEMENTS = 4096*4096 = 16777216
	sparsity 50%
	partition size 8
*/
int totalElements = 128*128;
int sparsity = 50;
int partSize = 8;


bool handleResultError(ostream& outfile, string name, BitMatrix& troubleMatrix, Matrix& resultTrue, Matrix& resultBad) {
	outfile<<endl;
	outfile<<name<<" code failed...\n";

	outfile<< "printing trouble matrix:\n";
	troubleMatrix.rprint(outfile);

	outfile<<"true answer:\n";
	resultTrue.rprint(outfile);

	outfile<<"wrong answer:\n";
	resultBad.rprint(outfile);
	return false;
}

void partitionTestf(int sparsity, const char* filename) {
	ofstream filez;
	cout<<"starting partition test\n";
	filez.open(filename);	

	filez << "sparsity 	partitions 	MSTSetupT 	MSTOps 	SEQSetupT 	SEQOps 	ORDOps\n";
	int columns = 128; 
	int rows = 128;
	int maxPartitions = 64;
	for(int i = 1; i < maxPartitions*2; i *= 2) {
		filez<<sparsity<<"	"<< i << "	";
		BitMatrix m(rows,columns); 
		Matrix v(columns,1);
		v.randomize(-1000,1000);
		m.setOneDensity(sparsity);

		Matrix resultMSTPart = m.multMSTEqPart(v, filez,i);
		Matrix resultSEQPart = m.multSequiterEqParts(v, filez, i);

		Matrix resultReg = m.mult(v, filez);

		assert(resultReg.eq(resultMSTPart) 
		|| handleResultError(filez, "multMSTPart", m, resultReg, resultMSTPart));

		assert(resultReg.eq(resultSEQPart) 
		|| handleResultError(filez, "multSequiterParts", m, resultReg, resultSEQPart));

		filez<<"\n";
	}

}

void sizeTestf(int sparsity, const char* filename) {
	ofstream filez;
	cout<<"starting size test\n";
	filez.open(filename);	

	filez << "sparsity 	size 	MSTSetupT 	MSTOps 	SEQSetupT 	SEQOps 	ORDOps\n";

	for(int i = 32; i < 1024; i *= 2) {
		int columns = i; 
		int rows = i;
		filez<<sparsity<<"	"<< i << "	";
		BitMatrix m(rows,columns); 
		Matrix v(columns,1);
		v.randomize(-1000,1000);
		m.setOneDensity(sparsity);

		Matrix resultMSTPart = m.multMSTEqPart(v, filez,16);
		Matrix resultSEQPart = m.multSequiterEqParts(v, filez, 1);

		Matrix resultReg = m.mult(v, filez);

		assert(resultReg.eq(resultMSTPart) 
		|| handleResultError(filez, "multMSTPart", m, resultReg, resultMSTPart));

		assert(resultReg.eq(resultSEQPart) 
		|| handleResultError(filez, "multSequiterParts", m, resultReg, resultSEQPart));

		filez<<"\n";
	}

}

void sparsityTestf(const char* filename) {
	ofstream filez;
	cout<<"starting size test\n";
	filez.open(filename);	

	filez << "sparsity 	MSTSetupT 	MSTOps 	SEQSetupT 	SEQOps 	ORDOps\n";
	int columns = 128; 
	int rows = 128;
	for(int i = 1; i < 96; i += 5) {
		
		filez<<i<< "	";
		BitMatrix m(rows,columns); 
		Matrix v(columns,1);
		v.randomize(-1000,1000);
		m.setOneDensity(i);

		Matrix resultMSTPart = m.multMSTEqPart(v, filez,16);
		Matrix resultSEQPart = m.multSequiterEqParts(v, filez, 1);

		Matrix resultReg = m.mult(v, filez);

		assert(resultReg.eq(resultMSTPart) 
		|| handleResultError(filez, "multMSTPart", m, resultReg, resultMSTPart));

		assert(resultReg.eq(resultSEQPart) 
		|| handleResultError(filez, "multSequiterParts", m, resultReg, resultSEQPart));

		filez<<"\n";
	}

}


void partitionTestVariance(int sparsity, char* filename, int parts) {
	ofstream filez;
	cout<<"starting partition test\n";
	filez.open(filename);
	filez << "sparsity 	partitions 	MST 	SEQ 	Ordinary\n";
	int columns = 128; 
	int rows = 128;
	int maxPartitions = 64;
	int i = parts;
	for(int j = 0; j < 25; j ++) {
		filez<<sparsity<<"	"<<i << "	";
		BitMatrix m(rows,columns); 
		Matrix v(columns,1);
		v.randomize(0,100);
		m.setOneDensity(sparsity);

		Matrix resultMSTPart = m.multMSTEqPart(v, filez,i);
		Matrix resultSEQPart = m.multSequiterEqParts(v, filez, i);

		Matrix resultReg = m.mult(v, filez);

		assert(resultReg.eq(resultMSTPart) 
		|| handleResultError(filez, "multMSTPart", m, resultReg, resultMSTPart));

		assert(resultReg.eq(resultSEQPart) 
		|| handleResultError(filez, "multSequiterParts", m, resultReg, resultSEQPart));
		filez<<"\n";
	}
}

void sizeTest() {
	ofstream filez;
	cout<<"starting size test\n";
	filez.open("sizeTest.txt");
	filez << "rows 	columns 	NumElements 	MST 	MSTpart 	SEQ 		SEQPART 	REG\n";
	for(int i = 0; i <5; i++) {
		cout<<"sizetest i = " << i << "/5" <<endl;
		int r=32<<i;
		int c=32<<i;
		filez << r << "	" << c << "	" << r*c << "	";
		BitMatrix m(r,c); 
		Matrix v(c,1);
		v.randomize(0,10);
		m.setOneDensity(50);

		Matrix resultMST = m.multMST(v, filez);
		Matrix resultMSTPart = m.multMSTPart(v, filez,8);
		Matrix resultSEQ = m.multSequiter(v, filez);
		Matrix resultSEQPart = m.multSequiterParts(v, filez, 8);

		Matrix resultReg = m.mult(v, filez);

		assert( resultReg.eq(resultMST) 
		|| handleResultError(filez, "multMST", m, resultReg, resultMST));

		assert(resultReg.eq(resultMSTPart) 
		|| handleResultError(filez, "multMSTPart", m, resultReg, resultMSTPart));

		assert(resultReg.eq(resultSEQ) 
		|| handleResultError(filez, "multSequiter", m, resultReg, resultSEQ));

		assert(resultReg.eq(resultSEQPart) 
		|| handleResultError(filez, "multSequiterParts", m, resultReg, resultSEQPart));

		filez<<"\n";
	}
	return;

}

void ratioTest2(int rows, int columns){
	ofstream filez;
	cout<<"Starting ratio test\n";
	filez.open("ratioTest.txt");
	filez << "rows 	columns 	NumElements 	MST 	MSTpart 	SEQ 		SEQPART 	REG\n";

}
void ratioTest() {
	ofstream filez;
	cout<<"Starting ratio test\n";
	filez.open("asd/ratioTest.txt");
	filez << "rows 	columns 	NumElements 	MST 	MSTpart 	SEQ 		SEQPART 	REG\n";
	int totalElements = 128*128;
	for(int i = 128; i <= totalElements/128; i *= 2) {
		cout<<"ratiotest i = " << i << "/16128" <<endl;
		int c = i;
		int r = totalElements/c;
		filez << r << "	" << c << "	" << r*c << "	";
		BitMatrix m(r,c); 
		Matrix v(c,1);
		v.randomize(0,10);
		m.setOneDensity(50);

		Matrix resultMSTPart = m.multMSTEqPart(v, filez,2);
		Matrix resultSEQPart = m.multSequiterEqParts(v, filez, 2);

		Matrix resultReg = m.mult(v, filez);

		assert(resultReg.eq(resultMSTPart) 
		|| handleResultError(filez, "multMSTPart", m, resultReg, resultMSTPart));

		assert(resultReg.eq(resultSEQPart) 
		|| handleResultError(filez, "multSequiterParts", m, resultReg, resultSEQPart));

		filez<<"\n";
	}
	return;

}

void sparsityTest() {
	ofstream filez;
	cout<<"starting sparsity test\n";
	filez.open("sparsityTest.txt");
	int r = 512;
	int c = 512;

	filez << "densityPct 	MST 	MSTpart 	SEQ 		SEQPART 	REG\n";
	for(int i = 1; i <= 100; i+=3) {
		filez << i << "	";
		cout<<"sparsityTest i = " << i << "/100" <<endl;
		BitMatrix m(r,c); 
		Matrix v(c,1);
		v.randomize(0,10);
		m.setOneDensity(i);

		Matrix resultMST = m.multMST(v, filez);
		Matrix resultMSTPart = m.multMSTPart(v, filez,8);
		Matrix resultSEQ = m.multSequiter(v, filez);
		Matrix resultSEQPart = m.multSequiterParts(v, filez, 8);

		Matrix resultReg = m.mult(v, filez);

		assert( resultReg.eq(resultMST) 
		|| handleResultError(filez, "multMST", m, resultReg, resultMST));

		assert(resultReg.eq(resultMSTPart) 
		|| handleResultError(filez, "multMSTPart", m, resultReg, resultMSTPart));

		assert(resultReg.eq(resultSEQ) 
		|| handleResultError(filez, "multSequiter", m, resultReg, resultSEQ));

		assert(resultReg.eq(resultSEQPart) 
		|| handleResultError(filez, "multSequiterParts", m, resultReg, resultSEQPart));

		filez<<"\n";
	}
	return;

}

void testParts(){
	for (int i=0; i<1000; i++)
	{
		int r = 64;
		int c = 64;
		BitMatrix m(r,c); 
		m.setOneDensity(50);
		Matrix v(c,1);
		v.randomize(0,10);

		m.buildSeqParts(8);
		m.buildSeq();

		Matrix resultSEQ = m.multSequiter(v, cout);
		Matrix resultSEQPart = m.multSequiterParts(v, cout, 8);
		assert(resultSEQPart.eq(resultSEQ) && "something is terribly wrong..");
	}
}


void MSTandSEQTiming(){
	// timing of MST, Seq, and the ordinary for a 128x128 matrix computing  
	// the 10000 power iterations for varying sparsity using the 
	// operation generation method


	ofstream timinglogfile;
	timinglogfile.open("timingsLog.dat");
	ofstream operandFileStream;
	ifstream MSTInStream; 
	ifstream SEQInStream;
	ifstream ORDInStream;

	string currentDirectory = "~/Research/S2017/";
	string SWPath = currentDirectory + "StopWatch.o";
	string operandFileName = "operand.dat";
	string operandFilePath = currentDirectory + operandFileName;

	string MSTFileName = "MSTfile.cpp"; 
	string MSTFilePath = currentDirectory + MSTFileName;
	string MSTExec = "MSTexec";
	string MSTExecPath = currentDirectory + MSTExec; 
	string MSTExecCommand = "" + MSTExecPath + " " + operandFilePath;
	string MSTOut = MSTExec+"_out";
	string MSTCompile = "g++ -w " +  MSTFilePath+ " " + SWPath + " -o "  + MSTExecPath; 


	string SEQFileName = "SEQfile.cpp"; 
	string SEQFilePath = currentDirectory + SEQFileName;
	string SEQExec = "SEQexec";
	string SEQExecPath = currentDirectory + SEQExec;
	string SEQExecCommand = "" + SEQExecPath + " " + operandFilePath; 
	string SEQOut = SEQExec+"_out";
	string SEQCompile = "g++ -w " +  SEQFilePath + " " + SWPath +" -o "  + SEQExecPath; 


	string ORDFileName = "ORDfile.cpp";
	string ORDFilePath = currentDirectory + ORDFileName;
	string ORDExec = "ORDexec";
	string ORDExecPath = currentDirectory + ORDExec; 
	string ORDExecCommand = "" + ORDExecPath + " " + operandFilePath;
	string ORDOut = ORDExec + "_out";
	string ORDCompile = "g++ -w " +  ORDFilePath + " " + SWPath +" -o "  + ORDExecPath; 


	int rows = 128; 
	int columns = 128; 



	for(int i = 1; i < 96; i +=5) {
		double MSTCycleCountSum = 0; 
		double SEQCycleCountSum = 0;
		double ORDCycleCountSum = 0; 
		vector<double> ORDtimings(25,0);
		vector<double> MSTtimings(25,0);
		vector<double> SEQtimings(25,0);
		for(int j = 0; j < 25; j ++) { // 25 tests for good STD
			BitMatrix mat(rows, columns); 
			Matrix operand(columns, 1);
			operand.randomize(-1000,1000);
			operandFileStream.open(operandFileName.c_str());
			for(int i=0; i<columns; i++) {
				operandFileStream << operand.getElement(i,0)<<endl;
			}
			operandFileStream.close();

			mat.setOneDensity(i);
			mat.buildORDProgram(ORDFileName.c_str());
			int err = 0; 
			err += system(ORDCompile.c_str());
			cout<<"compile ORD: " << err <<endl;
			err += system(ORDExecCommand.c_str());
			cout<<"run ORD: "<<err <<endl;	

			mat.buildMSTProgram(MSTFileName.c_str(), 16);
			// err += system(MSTCompile.c_str());
			cout<<"compile MST: " << err <<endl;
			err += system(MSTExecCommand.c_str());
			cout<<"run MST: "<< err <<endl;	

			mat.buildSEQProgram(SEQFileName.c_str());
			err += system(SEQCompile.c_str());
			cout<<"compile SEQ: " << err <<endl;
			err += system(SEQExecCommand.c_str());
			cout<<"run SEQ: "<<err <<endl;	
			if(err!=0){
				return; 
			}
			
			ORDInStream.open(ORDOut.c_str());
			double ORDCycles; 
			Matrix ORDResult(columns, 1);

			ORDInStream >> ORDCycles; 
			ORDCycleCountSum += ORDCycles;
			ORDtimings[j] = ORDCycles;

			for(int k = 0; k < columns; k++) {
				int element;
				ORDInStream >> element; 
				ORDResult.set(k,0, element);
			}
			ORDInStream.close();

			MSTInStream.open(MSTOut.c_str());
			double MSTCycles; 
			Matrix MSTResult(columns, 1);
			
			MSTInStream >> MSTCycles; 
			MSTCycleCountSum += MSTCycles;
			MSTtimings[j] = MSTCycles;


			for(int k = 0; k < columns; k++) {
				int element;
				MSTInStream >> element; 
				MSTResult.set(k,0, element);
			}
			MSTInStream.close();
			
			SEQInStream.open(SEQOut.c_str());
			double SEQCycles; 
			Matrix SEQResult(columns, 1);
			
			SEQInStream >> SEQCycles; 
			SEQCycleCountSum += SEQCycles;
			SEQtimings[j] = SEQCycles;


			for(int k = 0; k < columns; k++) {
				int element;
				SEQInStream >> element; 
				SEQResult.set(k,0, element);
			}
			SEQInStream.close();
			cout<<"total time ORD: " << ORDCycles << endl;
			cout<<"total time MST: " << MSTCycles << endl;
			cout<<"total time SEQ: " << SEQCycles << endl;
			cout<<"total time ORDsum: " << ORDCycleCountSum << endl;
			cout<<"total time MSTsum: " << MSTCycleCountSum << endl;
			cout<<"total time SEQsum: " << SEQCycleCountSum << endl;
			assert(MSTResult.eq(SEQResult) && "the results should be the same..");
		}
		double aveORD = ORDCycleCountSum/25;
		double aveMST = MSTCycleCountSum/25;
		double aveSEQ = SEQCycleCountSum/25;
		double ORDvar = 0;
		double MSTvar = 0;
		double SEQvar = 0;

		for (int j=0; j<25; j++) {
			double diff = ORDtimings[j] - aveORD;
			ORDvar += diff * diff;
			diff = ORDtimings[j] - aveORD;
			ORDvar += diff * diff;
			diff = MSTtimings[j] - aveMST;
			MSTvar += diff * diff;
			diff = SEQtimings[j] - aveSEQ;
			SEQvar += diff * diff;
		}
		ORDvar = ORDvar/25;
		MSTvar = MSTvar/25;
		SEQvar = SEQvar/25;

		double ORDstd = sqrt(ORDvar);
		double MSTstd = sqrt(MSTvar);
		double SEQstd = sqrt(SEQvar);


		cout<<"average time ORD: " << ORDCycleCountSum << endl;
		cout<<"average time MST: " << MSTCycleCountSum << endl;
		cout<<"average time SEQ: " << SEQCycleCountSum << endl;
		timinglogfile 
		<< i <<" " 
		<< aveORD << " " << ORDstd << " " 
		<< aveMST << " " << MSTstd << " " 
		<< aveSEQ << " " << SEQstd << "\n";
	}	
}


void opCounters(){
	for(int i=0; i<25; i++) {
		cout<<"*****************************************\n";
		cout<<"####### PASS "<<i+1<<"/25 ##########\n";
		cout<<"*****************************************\n";
		string fileLocation;
		string filePath;

		fileLocation = "data/sparsityTest/_";
		filePath = fileLocation + "sparsity_"+to_string(i)+".dat";
		sparsityTestf(filePath.c_str());

		fileLocation = "data/sizeTest/";
		filePath = fileLocation + "size.10_"+to_string(i)+".dat";
		sizeTestf(10, filePath.c_str());
		filePath = fileLocation + "size.25_"+to_string(i)+".dat";
		sizeTestf(25, filePath.c_str());
		filePath = fileLocation + "size.50_"+to_string(i)+".dat";
		sizeTestf(50, filePath.c_str());
		filePath = fileLocation + "size.75_"+to_string(i)+".dat";
		sizeTestf(75, filePath.c_str());
		filePath = fileLocation + "size.90_"+to_string(i)+".dat";
		sizeTestf(90, filePath.c_str());

		fileLocation = "data/partitionTest/";
		filePath = fileLocation + "part.10_"+to_string(i)+".dat";
		partitionTestf(10, filePath.c_str());
		filePath = fileLocation + "part.25_"+to_string(i)+".dat";
		partitionTestf(25, filePath.c_str());
		filePath = fileLocation + "part.50_"+to_string(i)+".dat";
		partitionTestf(50, filePath.c_str());
		filePath = fileLocation + "part.75_"+to_string(i)+".dat";
		partitionTestf(75, filePath.c_str());
		filePath = fileLocation + "part.90_"+to_string(i)+".dat";
	 	partitionTestf(90, filePath.c_str());
	}
}

void calcAvStdi(vector<int> vec, ostream& outfile) {
	int sum =0; 
	for(int i=0; i< vec.size(); i++){
		sum+=vec[i];
	}
	double average = sum / vec.size();
	double var = 0;
	for(int i=0; i< vec.size(); i++) {
		double diff = vec[i] - average;
		var += diff * diff;
	}
	var = var / vec.size();

	double std = sqrt(var);

	outfile << average << "	" << std << "	";
}

void calcAvStdd(vector<double> vec, ostream& outfile) {
	double sum =0; 
	for(int i=0; i< vec.size(); i++){
		sum+=vec[i]/(double)CLOCKS_PER_SEC;
	}
	double average = sum / vec.size();
	double var = 0;
	for(int i=0; i< vec.size(); i++) {
		double diff = vec[i]/(double)CLOCKS_PER_SEC - average;
		var += diff * diff;
	}
	var = var / vec.size();

	double std = sqrt(var);

	outfile << average << "	" << std << "	";
}

void opCountAvStd(){
	vector<ifstream> files(25);
	string fileLocation;
	string filePath;
	vector<double> MSTBuildT(25);
	vector<int> MSTOps(25);
	vector<double> SEQBuildT(25);
	vector<int> SEQOps(25);
	vector<int> ORDOps(25);

	ofstream sparsityTestFinal;
	sparsityTestFinal.open("data/sparsityTest/sparsityFinal");
	fileLocation = "data/sparsityTest/_";
	string dummy;
	for(int i=0; i<25; i++) {
		filePath = fileLocation + "sparsity_"+to_string(i)+".dat";
		files[i].open(filePath.c_str());
		for(int j=0; j<6; j++) {
			files[i] >> dummy;
		}
	}

	for(int j=0; j<19; j++) {
		int sparsity;
		for(int i=0; i<25; i++) {
			files[i] >> sparsity;
			files[i] >> MSTBuildT[i];
			files[i] >> MSTOps[i];
			files[i] >> SEQBuildT[i];
			files[i] >> SEQOps[i];
			files[i] >> ORDOps[i];
		}

		sparsityTestFinal << sparsity << "	";
		calcAvStdd(MSTBuildT, sparsityTestFinal);
		calcAvStdi(MSTOps, sparsityTestFinal);
		calcAvStdd(SEQBuildT, sparsityTestFinal);
		calcAvStdi(SEQOps, sparsityTestFinal);
		calcAvStdi(ORDOps, sparsityTestFinal);
		sparsityTestFinal<<endl;
	}

	for(int i=0; i<25; i++) {
		files[i].close();
	}
	
	sparsityTestFinal.close();
	//////////////////////////////////////////////
	//////////////////////////////////////////////
	ofstream partitionTestFinal10;
	partitionTestFinal10.open("data/partitionTest/partitionFinal10");
	fileLocation = "data/partitionTest/";
	
	for(int i=0; i<25; i++) {
		filePath = fileLocation + "part.10_"+to_string(i)+".dat";
		files[i].open(filePath.c_str());
		for(int j=0; j<7; j++) {
			files[i] >> dummy;
		}
	}

	for(int j=0; j<7; j++) {
		int sparsity;
		int partitions;
		for(int i=0; i<25; i++) {
			files[i] >> sparsity;
			files[i] >> partitions;
			files[i] >> MSTBuildT[i];
			files[i] >> MSTOps[i];
			files[i] >> SEQBuildT[i];
			files[i] >> SEQOps[i];
			files[i] >> ORDOps[i];
		}

		partitionTestFinal10 << partitions << "	";
		calcAvStdd(MSTBuildT, partitionTestFinal10);
		calcAvStdi(MSTOps, partitionTestFinal10);
		calcAvStdd(SEQBuildT, partitionTestFinal10);
		calcAvStdi(SEQOps, partitionTestFinal10);
		calcAvStdi(ORDOps, partitionTestFinal10);
		partitionTestFinal10<<endl;
	}
	
	for(int i=0; i<25; i++) {
		files[i].close();
	}
	partitionTestFinal10.close();
	//////////////////////////////////////
	ofstream partitionTestFinal25;
	partitionTestFinal25.open("data/partitionTest/partitionFinal25");
	fileLocation = "data/partitionTest/";
	
	for(int i=0; i<25; i++) {
		filePath = fileLocation + "part.25_"+to_string(i)+".dat";
		files[i].open(filePath.c_str());
		for(int j=0; j<7; j++) {
			files[i] >> dummy;
		}
	}

	for(int j=0; j<7; j++) {
		int sparsity;
		int partitions;
		for(int i=0; i<25; i++) {
			files[i] >> sparsity;
			files[i] >> partitions;
			files[i] >> MSTBuildT[i];
			files[i] >> MSTOps[i];
			files[i] >> SEQBuildT[i];
			files[i] >> SEQOps[i];
			files[i] >> ORDOps[i];
		}

		partitionTestFinal25 << partitions << "	";
		calcAvStdd(MSTBuildT, partitionTestFinal25);
		calcAvStdi(MSTOps, partitionTestFinal25);
		calcAvStdd(SEQBuildT, partitionTestFinal25);
		calcAvStdi(SEQOps, partitionTestFinal25);
		calcAvStdi(ORDOps, partitionTestFinal25);
		partitionTestFinal25<<endl;
	}
	
	for(int i=0; i<25; i++) {
		files[i].close();
	}
	partitionTestFinal25.close();
	////////////////////////////////
	ofstream partitionTestFinal50;
	partitionTestFinal50.open("data/partitionTest/partitionFinal50");
	fileLocation = "data/partitionTest/";
	
	for(int i=0; i<25; i++) {
		filePath = fileLocation + "part.50_"+to_string(i)+".dat";
		files[i].open(filePath.c_str());
		for(int j=0; j<7; j++) {
			files[i] >> dummy;
		}
	}

	for(int j=0; j<7; j++) {
		int sparsity;
		int partitions;
		for(int i=0; i<25; i++) {
			files[i] >> sparsity;
			files[i] >> partitions;
			files[i] >> MSTBuildT[i];
			files[i] >> MSTOps[i];
			files[i] >> SEQBuildT[i];
			files[i] >> SEQOps[i];
			files[i] >> ORDOps[i];
		}

		partitionTestFinal50 << partitions << "	";
		calcAvStdd(MSTBuildT, partitionTestFinal50);
		calcAvStdi(MSTOps, partitionTestFinal50);
		calcAvStdd(SEQBuildT, partitionTestFinal50);
		calcAvStdi(SEQOps, partitionTestFinal50);
		calcAvStdi(ORDOps, partitionTestFinal50);
		partitionTestFinal50<<endl;
	}
	
	for(int i=0; i<25; i++) {
		files[i].close();
	}
	partitionTestFinal50.close();
	///////////////////////////////////
	ofstream partitionTestFinal75;
	partitionTestFinal75.open("data/partitionTest/partitionFinal75");
	fileLocation = "data/partitionTest/";
	
	for(int i=0; i<25; i++) {
		filePath = fileLocation + "part.75_"+to_string(i)+".dat";
		files[i].open(filePath.c_str());
		for(int j=0; j<7; j++) {
			files[i] >> dummy;
		}
	}

	for(int j=0; j<7; j++) {
		int sparsity;
		int partitions;
		for(int i=0; i<25; i++) {
			files[i] >> sparsity;
			files[i] >> partitions;
			files[i] >> MSTBuildT[i];
			files[i] >> MSTOps[i];
			files[i] >> SEQBuildT[i];
			files[i] >> SEQOps[i];
			files[i] >> ORDOps[i];
		}

		partitionTestFinal75 << partitions << "	";
		calcAvStdd(MSTBuildT, partitionTestFinal75);
		calcAvStdi(MSTOps, partitionTestFinal75);
		calcAvStdd(SEQBuildT, partitionTestFinal75);
		calcAvStdi(SEQOps, partitionTestFinal75);
		calcAvStdi(ORDOps, partitionTestFinal75);
		partitionTestFinal75<<endl;
	}
	
	for(int i=0; i<25; i++) {
		files[i].close();
	}
	partitionTestFinal75.close();
	///////////////////////////////////////
	ofstream partitionTestFinal90;
	partitionTestFinal90.open("data/partitionTest/partitionFinal90");
	fileLocation = "data/partitionTest/";
	
	for(int i=0; i<25; i++) {
		filePath = fileLocation + "part.90_"+to_string(i)+".dat";
		files[i].open(filePath.c_str());
		for(int j=0; j<7; j++) {
			files[i] >> dummy;
		}
	}

	for(int j=0; j<7; j++) {
		int sparsity;
		int partitions;
		for(int i=0; i<25; i++) {
			files[i] >> sparsity;
			files[i] >> partitions;
			files[i] >> MSTBuildT[i];
			files[i] >> MSTOps[i];
			files[i] >> SEQBuildT[i];
			files[i] >> SEQOps[i];
			files[i] >> ORDOps[i];
		}

		partitionTestFinal90<< partitions << "	";
		calcAvStdd(MSTBuildT, partitionTestFinal90);
		calcAvStdi(MSTOps, partitionTestFinal90);
		calcAvStdd(SEQBuildT, partitionTestFinal90);
		calcAvStdi(SEQOps, partitionTestFinal90);
		calcAvStdi(ORDOps, partitionTestFinal90);
		partitionTestFinal90<<endl;
	}
	
	for(int i=0; i<25; i++) {
		files[i].close();
	}
	partitionTestFinal90.close();
	////////////////////////////////////////////////////
	///////////////////////////////////////////////////
	ofstream sizeTestFinal10;
	sizeTestFinal10.open("data/sizeTest/sizeFinal10");
	fileLocation = "data/sizeTest/";
	
	for(int i=0; i<25; i++) {
		filePath = fileLocation + "size.10_"+to_string(i)+".dat";
		files[i].open(filePath.c_str());
		for(int j=0; j<7; j++) {
			files[i] >> dummy;
		}
	}

	for(int j=0; j<5; j++) {
		int sparsity;
		int size;
		for(int i=0; i<25; i++) {
			files[i] >> sparsity;
			files[i] >> size;
			files[i] >> MSTBuildT[i];
			files[i] >> MSTOps[i];
			files[i] >> SEQBuildT[i];
			files[i] >> SEQOps[i];
			files[i] >> ORDOps[i];
		}

		sizeTestFinal10 << size << "	";
		calcAvStdd(MSTBuildT, sizeTestFinal10);
		calcAvStdi(MSTOps, sizeTestFinal10);
		calcAvStdd(SEQBuildT, sizeTestFinal10);
		calcAvStdi(SEQOps, sizeTestFinal10);
		calcAvStdi(ORDOps, sizeTestFinal10);
		sizeTestFinal10<<endl;
	}
	
	for(int i=0; i<25; i++) {
		files[i].close();
	}
	sizeTestFinal10.close();
	//////////////////////////////////////
	ofstream sizeTestFinal25;
	sizeTestFinal25.open("data/sizeTest/sizeFinal25");
	fileLocation = "data/sizeTest/";
	
	for(int i=0; i<25; i++) {
		filePath = fileLocation + "size.25_"+to_string(i)+".dat";
		files[i].open(filePath.c_str());
		for(int j=0; j<7; j++) {
			files[i] >> dummy;
		}
	}

	for(int j=0; j<5; j++) {
		int sparsity;
		int size;
		for(int i=0; i<25; i++) {
			files[i] >> sparsity;
			files[i] >> size;
			files[i] >> MSTBuildT[i];
			files[i] >> MSTOps[i];
			files[i] >> SEQBuildT[i];
			files[i] >> SEQOps[i];
			files[i] >> ORDOps[i];
		}

		sizeTestFinal25 << size << "	";
		calcAvStdd(MSTBuildT, sizeTestFinal25);
		calcAvStdi(MSTOps, sizeTestFinal25);
		calcAvStdd(SEQBuildT, sizeTestFinal25);
		calcAvStdi(SEQOps, sizeTestFinal25);
		calcAvStdi(ORDOps, sizeTestFinal25);
		sizeTestFinal25<<endl;
	}
	
	for(int i=0; i<25; i++) {
		files[i].close();
	}
	sizeTestFinal25.close();
	//////////////////////////////////////
	ofstream sizeTestFinal50;
	sizeTestFinal50.open("data/sizeTest/sizeFinal50");
	fileLocation = "data/sizeTest/";
	
	for(int i=0; i<25; i++) {
		filePath = fileLocation + "size.50_"+to_string(i)+".dat";
		files[i].open(filePath.c_str());
		for(int j=0; j<7; j++) {
			files[i] >> dummy;
		}
	}

	for(int j=0; j<5; j++) {
		int sparsity;
		int size;
		for(int i=0; i<25; i++) {
			files[i] >> sparsity;
			files[i] >> size;
			files[i] >> MSTBuildT[i];
			files[i] >> MSTOps[i];
			files[i] >> SEQBuildT[i];
			files[i] >> SEQOps[i];
			files[i] >> ORDOps[i];
		}

		sizeTestFinal50 << size << "	";
		calcAvStdd(MSTBuildT, sizeTestFinal50);
		calcAvStdi(MSTOps, sizeTestFinal50);
		calcAvStdd(SEQBuildT, sizeTestFinal50);
		calcAvStdi(SEQOps, sizeTestFinal50);
		calcAvStdi(ORDOps, sizeTestFinal50);
		sizeTestFinal50<<endl;
	}
	
	for(int i=0; i<25; i++) {
		files[i].close();
	}
	sizeTestFinal50.close();
	//////////////////////////////////////
	ofstream sizeTestFinal75;
	sizeTestFinal75.open("data/sizeTest/sizeFinal75");
	fileLocation = "data/sizeTest/";
	
	for(int i=0; i<25; i++) {
		filePath = fileLocation + "size.75_"+to_string(i)+".dat";
		files[i].open(filePath.c_str());
		for(int j=0; j<7; j++) {
			files[i] >> dummy;
		}
	}

	for(int j=0; j<5; j++) {
		int sparsity;
		int size;
		for(int i=0; i<25; i++) {
			files[i] >> sparsity;
			files[i] >> size;
			files[i] >> MSTBuildT[i];
			files[i] >> MSTOps[i];
			files[i] >> SEQBuildT[i];
			files[i] >> SEQOps[i];
			files[i] >> ORDOps[i];
		}

		sizeTestFinal75 << size << "	";
		calcAvStdd(MSTBuildT, sizeTestFinal75);
		calcAvStdi(MSTOps, sizeTestFinal75);
		calcAvStdd(SEQBuildT, sizeTestFinal75);
		calcAvStdi(SEQOps, sizeTestFinal75);
		calcAvStdi(ORDOps, sizeTestFinal75);
		sizeTestFinal75<<endl;
	}
	
	for(int i=0; i<25; i++) {
		files[i].close();
	}
	sizeTestFinal75.close();
	//////////////////////////////////////
	ofstream sizeTestFinal90;
	sizeTestFinal90.open("data/sizeTest/sizeFinal90");
	fileLocation = "data/sizeTest/";
	
	for(int i=0; i<25; i++) {
		filePath = fileLocation + "size.90_"+to_string(i)+".dat";
		files[i].open(filePath.c_str());
		for(int j=0; j<7; j++) {
			files[i] >> dummy;
		}
	}

	for(int j=0; j<5; j++) {
		int sparsity;
		int size;
		for(int i=0; i<25; i++) {
			files[i] >> sparsity;
			files[i] >> size;
			files[i] >> MSTBuildT[i];
			files[i] >> MSTOps[i];
			files[i] >> SEQBuildT[i];
			files[i] >> SEQOps[i];
			files[i] >> ORDOps[i];
		}

		sizeTestFinal90 << size << "	";
		calcAvStdd(MSTBuildT, sizeTestFinal90);
		calcAvStdi(MSTOps, sizeTestFinal90);
		calcAvStdd(SEQBuildT, sizeTestFinal90);
		calcAvStdi(SEQOps, sizeTestFinal90);
		calcAvStdi(ORDOps, sizeTestFinal90);
		sizeTestFinal90<<endl;
	}
	
	for(int i=0; i<25; i++) {
		files[i].close();
	}
	sizeTestFinal90.close();
	//////////////////////////////////////
	

}

int main(int argc, char *argv[])
{	
	if(argc==1) {


			
			//partitionTest(50, "part.timeBigO3.50.3-1.txt");
			//ratioTest();
			//testFileWriting();
			//testSeqFile();
			//MSTandSEQTiming();
			opCountAvStd();
			return 0;

	}

	else {
		string input; 
		enum options {ex, size, sparsity, ratio, all};
		std::unordered_map<string, unsigned int> optionsMap = {
			{"exit", ex},
			{"size", size},
			{"sparsity", sparsity},
			{"ratio", ratio},
			{"all", all}
		};


		bool running = true;

		while(running) {
			cin>>input;
			switch(optionsMap[input]) {
				case all: break; 
				case size: sizeTest(); break;
				case sparsity: sparsityTest(); break;
				case ratio: ratioTest(); break;
				case ex: running = false; break;
			}
		}
		return 0;
	}
	
}