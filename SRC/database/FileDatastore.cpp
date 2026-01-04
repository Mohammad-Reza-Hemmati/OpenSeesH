/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.20 $
// $Date: 2007-04-05 01:15:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/FileDatastore.cpp,v $


// Written: fmk 
// Created: 10/98
//
// Description: This file contains the class implementation for FileDatastore.
// FileDatastore is a concrete subclas of FE_Datastore. A FileDatastore 
// object is used in the program to store/restore the geometry and state 
// information in a domain at a particular instance in the analysis. The
// information is stored in text files.
//
// What: "@(#) FileDatastore.C, revA"

#include "FileDatastore.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <bool.h>
#include <iostream>
#include <iomanip>
using std::ios;
using std::setiosflags;


#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include <ostream>   // for std::ostream
#include <set>

#define _DEBUG

#ifdef _DEBUG

/* PSEUDOCODE / DETAILED PLAN
 - Maintain a single shared std::ofstream per dataBase (keyed by dataBase string).
 - Maintain a per-dataBase set<int> tracking which dataTags have already had a header written.
 - getDataTagStream(dataTag, dataBase):
		 * Compute key from dataBase (empty string if null).
		 * Look up existing ofstream* in dataBaseStreams.
		 * If not present:
				 - Try to open "<dataBase>/debug/datastore.log" in append mode.
				 - If that fails, fall back to "<dataBase>.datastore.log".
				 - If both fail, return nullptr.
				 - Store the opened stream in dataBaseStreams[key].
		 * Using dataBaseHeaderWritten[key] (auto-created if absent) check if dataTag is present.
				 - If dataTag not present, write a single header line "=== dataTag: <n> ===\n" and flush, then mark dataTag written.
				 - If dataTag already present, do not write another header.
		 * Return the ofstream* (may be used many times for the same dataTag; header only on first write).
 - closeStreamsForDataBase(dataBase):
		 * If dataBase is null return.
		 * Lookup stream in dataBaseStreams and, if present, close and delete the stream and erase map entry.
		 * Erase corresponding header-tracking set in dataBaseHeaderWritten to free memory.
 - Keep helper functions for pretty-printing ID/Vector/Matrix (writeIDToStream, writeVectorToStream, writeMatrixToStream).
 - This ensures multiple writes for the same dataTag use the same file and only the first write is preceded by a headline.
*/

namespace {

	// buffer per dataBase and per dataTag; map keeps keys and tags sorted (ascending)
	static std::map<std::string, std::map<int, std::vector<std::string>>> dataBaseBuffers;

	// Lazily open stream for key; tries <key>/debug/datastore.log then <key>.datastore.log
	static std::ofstream* openStreamForKey(const std::string& key) {
		std::ostringstream fname;
		fname << key << "/debug/datastore.log";
		std::ofstream* ofs = new std::ofstream(fname.str().c_str(), std::ios::out | std::ios::app);
		return ofs;
	}

	// Append a formatted line to the in-memory buffer for dataTag/key.
	static void appendLogEntry(int dataTag, const char* dataBase, const std::string& line) {
		dataBaseBuffers[dataBase][dataTag].push_back(line);
	}

	// Flush buffered entries for a dataBase to its log file in ascending order of dataTag.
	static void flushBuffersForDataBase(const char* dataBase) {
		auto bufIt = dataBaseBuffers.find(dataBase);
		if (bufIt == dataBaseBuffers.end())
			return; // nothing buffered

		std::ofstream* ofs = openStreamForKey(dataBase);
		if (ofs == nullptr)
			return; // cannot open log file; drop entries

		auto& tagMap = bufIt->second; // std::map<int, vector<string>> ensures ascending order
		for (auto const& tagPair : tagMap) {
			int dataTag = tagPair.first;
			const std::vector<std::string>& lines = tagPair.second;
			(*ofs) << "=== dataTag: " << dataTag << " ===\n";
			for (const auto& l : lines) {
				(*ofs) << l << "\n";
			}
		}
		ofs->flush();
		// clear buffer for this key
		dataBaseBuffers.erase(bufIt);
	}

	// Existing helper formatting functions write to any ostream (ostringstream supported)
	static void writeIDToStream(std::ostream& os, const ID& id) {
		int sz = id.Size();
		os << "[";
		for (int i = 0; i < sz; ++i) {
			if (i) os << " ";
			os << id(i);
		}
		os << "]";
	}

	static void writeVectorToStream(std::ostream& os, const Vector& v) {
		int sz = v.Size();
		os << "[";
		for (int i = 0; i < sz; ++i) {
			if (i) os << " ";
			os << std::setprecision(16) << v(i);
		}
		os << "]";
	}

	static void writeMatrixToStream(std::ostream& os, const Matrix& m) {
		int rows = m.noRows();
		int cols = m.noCols();
		os << "[";
		for (int r = 0; r < rows; ++r) {
			if (r) os << " ";
			os << "[";
			for (int c = 0; c < cols; ++c) {
				if (c) os << " ";
				os << std::setprecision(16) << m(r, c);
			}
			os << "]";
		}
		os << "]";
	}

}
#endif

FileDatastore::FileDatastore(const char* dataBaseName,
	Domain& theDomain,
	FEM_ObjectBroker& theObjBroker)
	:FE_Datastore(theDomain, theObjBroker),
	data(0), sizeData(0), currentMaxInt(0), currentMaxDouble(0), currentCommitTag(-1)
{
	dataBase = new char[strlen(dataBaseName) + 1];
	strcpy(dataBase, dataBaseName);

	this->resizeDouble(1024);
}

void FileDatastore::Flush()
{
#ifdef _DEBUG
	flushBuffersForDataBase(dataBase);
#endif
}

FileDatastore::~FileDatastore()
{
	if (dataBase != 0)
		delete[] dataBase;

	if (data != 0)
		delete[] data;


	//  while (theIDFilesIter != theIDFiles.end()) {
	//    theIDFilesIter++;
	// }

	for (theIDFilesIter = theIDFiles.begin(); theIDFilesIter != theIDFiles.end(); theIDFilesIter++) {
		FileDatastoreOutputFile* theFileStruct = theIDFilesIter->second;
		fstream* theFile = theFileStruct->theFile;
		if (theFile != 0) {
			theFile->close();
			delete theFile;
		}
		delete theFileStruct;
	}

	theIDFiles.clear();

	for (theMatFilesIter = theMatFiles.begin(); theMatFilesIter != theMatFiles.end(); theMatFilesIter++) {
		FileDatastoreOutputFile* theFileStruct = theMatFilesIter->second;
		fstream* theFile = theFileStruct->theFile;
		if (theFile != 0) {
			theFile->close();
			delete theFile;
		}
		delete theFileStruct;
	}
	theMatFiles.clear();
	for (theVectFilesIter = theVectFiles.begin(); theVectFilesIter != theVectFiles.end(); theVectFilesIter++) {
		FileDatastoreOutputFile* theFileStruct = theVectFilesIter->second;
		fstream* theFile = theFileStruct->theFile;
		if (theFile != 0) {
			theFile->close();
			delete theFile;
		}
		delete theFileStruct;
	}
	theVectFiles.clear();
}


int
FileDatastore::commitState(int commitTag)
{
	int result = FE_Datastore::commitState(commitTag);
	if (result == commitTag)
		resetFilePointers();

	return result;
}


void
FileDatastore::resetFilePointers(void) {
	for (theIDFilesIter = theIDFiles.begin(); theIDFilesIter != theIDFiles.end(); theIDFilesIter++) {
		FileDatastoreOutputFile* theFileStruct = theIDFilesIter->second;
		fstream* theFile = theFileStruct->theFile;
		if (theFile != 0) {
			theFile->seekp(0, ios::beg);
			*(theIntData.dbTag) = theFileStruct->maxDbTag;
			theFile->write(data, sizeof(int));
			theFile->close();
			delete theFile;
			theFileStruct->theFile = 0;
		}
	}

	for (theMatFilesIter = theMatFiles.begin(); theMatFilesIter != theMatFiles.end(); theMatFilesIter++) {
		FileDatastoreOutputFile* theFileStruct = theMatFilesIter->second;
		fstream* theFile = theFileStruct->theFile;
		if (theFile != 0) {
			theFile->seekp(0, ios::beg);
			*(theIntData.dbTag) = theFileStruct->maxDbTag;
			theFile->write(data, sizeof(int));
			theFile->close();
			delete theFile;
			theFileStruct->theFile = 0;
		}
	}

	for (theVectFilesIter = theVectFiles.begin(); theVectFilesIter != theVectFiles.end(); theVectFilesIter++) {
		FileDatastoreOutputFile* theFileStruct = theVectFilesIter->second;
		fstream* theFile = theFileStruct->theFile;
		if (theFile != 0) {
			theFile->seekp(0, ios::beg);
			*(theIntData.dbTag) = theFileStruct->maxDbTag;
			theFile->write(data, sizeof(int));
			theFile->close();
			delete theFile;
			theFileStruct->theFile = 0;
		}
	}
	currentCommitTag = -1;
}


int
FileDatastore::sendMsg(int dataTag, int commitTag,
	const Message&,
	ChannelAddress* theAddress)
{
	opserr << "FileDatastore::sendMsg() - not yet implemented\n";
	return -1;
}

int
FileDatastore::recvMsg(int dataTag, int commitTag,
	Message&,
	ChannelAddress* theAddress)
{
	opserr << "FileDatastore::recvMsg() - not yet implemented\n";
	return -1;
}

int
FileDatastore::recvMsgUnknownSize(int dataTag, int commitTag,
	Message&,
	ChannelAddress* theAddress)
{
	opserr << "FileDatastore::recvMsgUnknownSize() - not yet implemented\n";
	return -1;
}


int
FileDatastore::sendID(int dataTag, int commitTag,
	const ID& theID,
	ChannelAddress* theAddress)
{
#ifdef _DEBUG
	// write to per-dataTag file instead of opserr
	std::ostringstream oss;
	oss << "ID: ";
	writeIDToStream(oss, theID);
	appendLogEntry(dataTag, dataBase, oss.str());
#endif
	if (currentCommitTag != commitTag)
		this->resetFilePointers();

	currentCommitTag = commitTag;

	FileDatastoreOutputFile* theFileStruct = 0;

	//
	// next we see if we already have this file; 
	//  if not we need to create data structure & open it
	//  if we have data structure, need to check file is opened (we close in a commit)
	//

	int idSize = theID.Size();
	int stepSize = (1 + idSize) * sizeof(int);

	theIDFilesIter = theIDFiles.find(idSize);
	if (theIDFilesIter == theIDFiles.end()) {

		// we first check if we need to resize send buffer
		if (idSize > currentMaxInt) {
			if (this->resizeInt(idSize) < 0) {
				opserr << "FileDatastore::sendID() - failed in resizeInt()\n";
				return -1;
			}
		}

		char* fileName = new char[strlen(dataBase) + 21];
		theFileStruct = new FileDatastoreOutputFile;

		if (fileName == 0 || theFileStruct == 0) {
			opserr << "FileDatastore::sendID() - out of memory\n";
			return -1;
		}

		static char intName[20];
		strcpy(fileName, dataBase);
		sprintf(intName, "%d.%d", idSize, commitTag);
		strcat(fileName, ".IDs.");
		strcat(fileName, intName);

		if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
			opserr << "FileDatastore::sendID() - could not open file\n";
			delete[] fileName;
			return -1;
		}
		else
			theIDFiles.insert(MAP_FILES_TYPE(idSize, theFileStruct));

		delete[] fileName;

	}
	else {

		theFileStruct = theIDFilesIter->second;

		// make sure not close from a last commit
		if (theFileStruct->theFile == 0) {

			if (idSize > currentMaxInt) {
				if (this->resizeInt(idSize) < 0) {
					opserr << "FileDatastore::sendID() - failed in resizeInt()\n";
					return -1;
				}
			}
			char* fileName = new char[strlen(dataBase) + 21];
			static char intName[20];
			strcpy(fileName, dataBase);
			sprintf(intName, "%d.%d", idSize, commitTag);
			strcat(fileName, ".IDs.");
			strcat(fileName, intName);

			if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
				opserr << "FileDatastore::sendID() - could not open file\n";
				delete[] fileName;
				return -1;

			}
			delete[] fileName;
		}
	}

	//
	// find location in file to place the data
	//

	fstream* theStream = theFileStruct->theFile;
	bool found = false;
	STREAM_POSITION_TYPE pos = theStream->tellg();
	STREAM_POSITION_TYPE fileEnd = theFileStruct->fileEnd;

	// we first check if the data can go at the end of the file
	// true if commitTag larger than any we have encountered so far

	if (theFileStruct->maxDbTag < dataTag) {
		pos = fileEnd;
		found = true;
		theFileStruct->maxDbTag = dataTag;
	}

	// try current location    
	if (pos < fileEnd && found == false) {
		// #ifdef _WIN32
				// must be a bug in the vc compiler! .. we are here already tellg() above!!
		theStream->seekg(pos, ios::beg);
		// #endif
		theStream->read(data, stepSize);
		if (*(theIntData.dbTag) == dataTag) {
			found = true;
		}
	}

	// we have to search from the beginning of the file  
	if (found == false) {
		*(theIntData.dbTag) = -1;
		pos = sizeof(int);
		theStream->seekg(pos, ios::beg);
		while (pos < fileEnd && found == false) {
			theStream->read(data, stepSize);

			if (*(theIntData.dbTag) == dataTag)
				found = true;
			else
				pos += stepSize;
		}

		if (found == false) {
			pos = fileEnd;
		}
	}

	//
	// we now place the data to be sent into our buffer
	//

	*(theIntData.dbTag) = dataTag;
	for (int i = 0; i < idSize; i++)
		theIntData.data[i] = theID(i);

	//
	// we now write the data
	//

	theStream->seekp(pos, ios::beg); // reset so can go write at the end

	theStream->write(data, stepSize);
	if (theStream->bad()) {
		opserr << "FileDatastore::sendID() - error writing to file\n";
		return -1;
	}

	// update the size of file if we have added to eof
	if (fileEnd <= pos)
		theFileStruct->fileEnd = pos + stepSize;

	//opserr << "WROTE: " << dataTag << " " << pos << endln;

	return 0;
}

int
FileDatastore::recvID(int dataTag, int commitTag,
	ID& theID,
	ChannelAddress* theAddress)
{
	if (currentCommitTag != commitTag)
		this->resetFilePointers();

	currentCommitTag = commitTag;

	FileDatastoreOutputFile* theFileStruct;

	//
	// next we see if we already have this file; 
	//  if not we need to create data structure & open it
	//  if we have data structure, need to check file is opened (we close in a commit)
	//

	int idSize = theID.Size();
	int stepSize = (1 + idSize) * sizeof(int);

	theIDFilesIter = theIDFiles.find(idSize);
	if (theIDFilesIter == theIDFiles.end()) {

		// we first check if we need to resize recv buffer
		if (idSize > currentMaxInt) {
			if (this->resizeInt(idSize) < 0) {
				opserr << "FileDatastore::recvID() - failed in resizeInt()\n";
				return -1;
			}
		}

		char* fileName = new char[strlen(dataBase) + 21];
		theFileStruct = new FileDatastoreOutputFile;

		if (fileName == 0 || theFileStruct == 0) {
			opserr << "FileDatastore::recvID() - out of memory\n";
			return -1;
		}

		static char intName[20];
		strcpy(fileName, dataBase);
		sprintf(intName, "%d.%d", idSize, commitTag);
		strcat(fileName, ".IDs.");
		strcat(fileName, intName);

		if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
			opserr << "FileDatastore::recvID() - could not open file\n";
			delete[] fileName;
			return -1;
		}
		else
			theIDFiles.insert(MAP_FILES_TYPE(idSize, theFileStruct));

		delete[] fileName;
	}
	else {

		theFileStruct = theIDFilesIter->second;

		// make sure not close from a last commit
		if (theFileStruct->theFile == 0) {

			if (idSize > currentMaxInt) {
				if (this->resizeInt(idSize) < 0) {
					opserr << "FileDatastore::recvID() - failed in resizeInt()\n";
					return -1;
				}
			}
			char* fileName = new char[strlen(dataBase) + 21];
			static char intName[20];
			strcpy(fileName, dataBase);
			sprintf(intName, "%d.%d", idSize, commitTag);
			strcat(fileName, ".IDs.");
			strcat(fileName, intName);

			if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
				opserr << "FileDatastore::recvID() - could not open file\n";
				delete[] fileName;
				return -1;
			}
			delete[] fileName;
		}
	}

	fstream* theStream = theFileStruct->theFile;
	STREAM_POSITION_TYPE fileEnd = theFileStruct->fileEnd;
	STREAM_POSITION_TYPE pos = theStream->tellg();

	//
	// find location in file to place the data
	//

	bool found = false;

	// we try the current file position first

	if (pos < fileEnd) {
		theStream->read(data, stepSize);
		if (*(theIntData.dbTag) == dataTag) {
			found = true;
			pos += stepSize;
		}
	}

	// we must search from the beginning of the file
	if (found == false) {
		pos = sizeof(int);
		theStream->seekg(pos, ios::beg);
		while ((pos < fileEnd) && (found == false)) {
			theStream->read(data, stepSize);
			if (*(theIntData.dbTag) == dataTag)
				found = true;
			else
				pos += stepSize;
		}
	}

	if (found == false) {
		opserr << "FileDatastore::recvID() - failed\n";
		return -1;
	}

	//opserr << "READ: " << dataTag << " " << pos << endln;

	// we now place the received data into the ID 
	for (int i = 0; i < idSize; i++)
		theID(i) = theIntData.data[i];

#ifdef _DEBUG
	std::ostringstream oss;
	oss << "ID: ";
	writeIDToStream(oss, theID);
	appendLogEntry(dataTag, dataBase, oss.str());
#endif
	return 0;
}

int
FileDatastore::sendMatrix(int dataTag, int commitTag,
	const Matrix& theMatrix,
	ChannelAddress* theAddress)
{
#ifdef _DEBUG
	std::ostringstream oss;
	oss << "Matrix: ";
	writeMatrixToStream(oss, theMatrix);
	appendLogEntry(dataTag, dataBase, oss.str());
#endif
	if (currentCommitTag != commitTag)
		this->resetFilePointers();

	currentCommitTag = commitTag;

	FileDatastoreOutputFile* theFileStruct;

	//
	// next we see if we already have this file; 
	//  if not we need to create data structure & open it
	//  if we have data structure, need to check file is opened (we close in a commit)
	//

	// we first ensure that the Matrix is not too big
	int noMatCols = theMatrix.noCols();
	int noMatRows = theMatrix.noRows();
	int matSize = noMatRows * noMatCols;;
	int stepSize = sizeof(int) + matSize * sizeof(double);

	theMatFilesIter = theMatFiles.find(matSize);
	if (theMatFilesIter == theMatFiles.end()) {

		// we first check if we need to resize send buffer
		if (matSize > currentMaxDouble) {
			if (this->resizeDouble(matSize) < 0) {
				opserr << "FileDatastore::sendMatrix() - failed in resizeInt()\n";
				return -1;
			}
		}

		char* fileName = new char[strlen(dataBase) + 21];
		theFileStruct = new FileDatastoreOutputFile;

		if (fileName == 0 || theFileStruct == 0) {
			opserr << "FileDatastore::sendMatrix() - out of memory\n";
			return -1;
		}

		static char intName[20];
		strcpy(fileName, dataBase);
		sprintf(intName, "%d.%d", matSize, commitTag);
		strcat(fileName, ".MATs.");
		strcat(fileName, intName);

		if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
			opserr << "FileDatastore::sendMatrix() - could not open file\n";
			delete[] fileName;
			return -1;
		}
		else
			theMatFiles.insert(MAP_FILES_TYPE(matSize, theFileStruct));

		delete[] fileName;

	}
	else {

		theFileStruct = theMatFilesIter->second;

		// make sure not close from a last commit
		if (theFileStruct->theFile == 0) {

			if (matSize > currentMaxDouble) {
				if (this->resizeDouble(matSize) < 0) {
					opserr << "FileDatastore::sendMatrix() - failed in resizeInt()\n";
					return -1;
				}
			}
			char* fileName = new char[strlen(dataBase) + 21];
			static char intName[20];
			strcpy(fileName, dataBase);
			sprintf(intName, "%d.%d", matSize, commitTag);
			strcat(fileName, ".MATs.");
			strcat(fileName, intName);

			if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
				opserr << "FileDatastore::sendMatrix() - could not open file\n";
				delete[] fileName;
				return -1;
			}

			delete[] fileName;
		}
	}

	//
	// find location in file to place the data
	//

	fstream* theStream = theFileStruct->theFile;
	bool found = false;
	STREAM_POSITION_TYPE pos = theStream->tellg();
	STREAM_POSITION_TYPE fileEnd = theFileStruct->fileEnd;

	// we first check if the data can go at the end of the file
	// true if commitTag larger than any we have encountered so far

	if (theFileStruct->maxDbTag < dataTag) {
		pos = fileEnd;
		found = true;
		theFileStruct->maxDbTag = dataTag;
	}

	// try current location    
	if (pos < fileEnd && found == false) {
		// #ifdef _WIN32
				// must be a bug in the vc compiler! .. we are here already tellg() above!!
		theStream->seekg(pos, ios::beg);
		// #endif
		theStream->read(data, stepSize);
		if (*(theIntData.dbTag) == dataTag) {
			found = true;
		}
	}

	// we have to search from the beginning of the file  
	if (found == false) {
		*(theIntData.dbTag) = -1;
		pos = sizeof(int);
		theStream->seekg(pos, ios::beg);
		while (pos < fileEnd && found == false) {
			theStream->read(data, stepSize);

			if (*(theIntData.dbTag) == dataTag)
				found = true;
			else
				pos += stepSize;
		}

		if (found == false) {
			pos = fileEnd;
		}
	}

	//
	// we now place the data to be sent into our buffer
	//

	*(theDoubleData.dbTag) = dataTag;
	int loc = 0;
	for (int j = 0; j < noMatCols; j++)
		for (int k = 0; k < noMatRows; k++) {
			theDoubleData.data[loc] = theMatrix(k, j);
			loc++;
		}

	//
	// we now write the data
	//

	theStream->seekp(pos, ios::beg); // reset so can go write at the end
	theStream->write(data, stepSize);

	// update the size of file if we have added to eof
	if (theFileStruct->fileEnd <= pos)
		theFileStruct->fileEnd += stepSize;

	return 0;
}

int
FileDatastore::recvMatrix(int dataTag, int commitTag,
	Matrix& theMatrix,
	ChannelAddress* theAddress)
{
	if (currentCommitTag != commitTag)
		this->resetFilePointers();

	currentCommitTag = commitTag;

	FileDatastoreOutputFile* theFileStruct;

	//
	// next we see if we already have this file; 
	//  if not we need to create data structure & open it
	//  if we have data structure, need to check file is opened (we close in a commit)
	//

	// we first check Matrix not too big
	int noMatCols = theMatrix.noCols();
	int noMatRows = theMatrix.noRows();
	int matSize = noMatRows * noMatCols;;
	int stepSize = sizeof(int) + matSize * sizeof(double);

	theMatFilesIter = theMatFiles.find(matSize);
	if (theMatFilesIter == theMatFiles.end()) {

		// we first check if we need to resize recv buffer
		if (matSize > currentMaxDouble) {
			if (this->resizeDouble(matSize) < 0) {
				opserr << "FileDatastore::recvMatrix() - failed in resizeDouble()\n";
				return -1;
			}
		}

		char* fileName = new char[strlen(dataBase) + 21];
		theFileStruct = new FileDatastoreOutputFile;

		if (fileName == 0 || theFileStruct == 0) {
			opserr << "FileDatastore::recvMatrix() - out of memory\n";
			return -1;
		}

		static char intName[20];
		strcpy(fileName, dataBase);
		sprintf(intName, "%d.%d", matSize, commitTag);
		strcat(fileName, ".MATs.");
		strcat(fileName, intName);

		if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
			opserr << "FileDatastore::recvMatrix() - could not open file\n";
			delete[] fileName;
			return -1;
		}
		else
			theMatFiles.insert(MAP_FILES_TYPE(matSize, theFileStruct));

		delete[] fileName;

	}
	else {

		theFileStruct = theMatFilesIter->second;

		// make sure not close from a last commit
		if (theFileStruct->theFile == 0) {

			if (matSize > currentMaxDouble) {
				if (this->resizeDouble(matSize) < 0) {
					opserr << "FileDatastore::recvMatrix() - failed in resizeInt()\n";
					return -1;
				}
			}
			char* fileName = new char[strlen(dataBase) + 21];
			static char intName[20];
			strcpy(fileName, dataBase);
			sprintf(intName, "%d.%d", matSize, commitTag);
			strcat(fileName, ".MATs.");
			strcat(fileName, intName);

			if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
				opserr << "FileDatastore::recvMatrix() - could not open file\n";
				delete[] fileName;
				return -1;
			}

			delete[] fileName;
		}
	}


	fstream* theStream = theFileStruct->theFile;
	STREAM_POSITION_TYPE fileEnd = theFileStruct->fileEnd;
	STREAM_POSITION_TYPE pos = theStream->tellg();

	//
	// find location in file to place the data
	//

	bool found = false;

	// we try the current file position first

	if (pos < fileEnd) {
		theStream->read(data, stepSize);
		if (*(theIntData.dbTag) == dataTag) {
			found = true;
			pos += stepSize;
		}
	}

	// we must search from the beginning of the file
	if (found == false) {
		pos = sizeof(int);
		theStream->seekg(pos, ios::beg);
		while ((pos < fileEnd) && (found == false)) {
			theStream->read(data, stepSize);
			if (*(theIntData.dbTag) == dataTag)
				found = true;
			else
				pos += stepSize;
		}
	}

	if (found == false) {
		opserr << "FileDatastore::recvMatrix() - failed\n";
		return -1;
	}

	int loc = 0;
	for (int j = 0; j < noMatCols; j++)
		for (int k = 0; k < noMatRows; k++) {
			theMatrix(k, j) = theDoubleData.data[loc];
			loc++;
		}

#ifdef _DEBUG
	std::ostringstream oss;
	oss << "Matrix: ";
	writeMatrixToStream(oss, theMatrix);
	appendLogEntry(dataTag, dataBase, oss.str());
#endif
	return 0;
}

int
FileDatastore::sendVector(int dataTag, int commitTag,
	const Vector& theVector,
	ChannelAddress* theAddress)
{
#ifdef _DEBUG
	std::ostringstream oss;
	oss << "Vector: ";
	writeVectorToStream(oss, theVector);
	appendLogEntry(dataTag, dataBase, oss.str());
#endif
	if (currentCommitTag != commitTag)
		this->resetFilePointers();

	currentCommitTag = commitTag;

	FileDatastoreOutputFile* theFileStruct;

	//
	// next we see if we already have this file; 
	//  if not we need to create data structure & open it
	//  if we have data structure, need to check file is opened (we close in a commit)
	//

	// we first ensure that the Matrix is not too big
	int vectSize = theVector.Size();
	int stepSize = sizeof(int) + vectSize * sizeof(double);

	theVectFilesIter = theVectFiles.find(vectSize);
	if (theVectFilesIter == theVectFiles.end()) {

		// we first check if we need to resize send buffer
		if (vectSize > currentMaxDouble) {
			if (this->resizeDouble(vectSize) < 0) {
				opserr << "FileDatastore::sendVector() - failed in resizeInt()\n";
				return -1;
			}
		}

		char* fileName = new char[strlen(dataBase) + 21];
		theFileStruct = new FileDatastoreOutputFile;

		if (fileName == 0 || theFileStruct == 0) {
			opserr << "FileDatastore::sendVector() - out of memory\n";
			return -1;
		}

		static char intName[20];
		strcpy(fileName, dataBase);
		sprintf(intName, "%d.%d", vectSize, commitTag);
		strcat(fileName, ".VECs.");
		strcat(fileName, intName);

		if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
			opserr << "FileDatastore::sendVector() - could not open file\n";
			delete[] fileName;
			return -1;
		}
		else
			theVectFiles.insert(MAP_FILES_TYPE(vectSize, theFileStruct));

		delete[] fileName;

	}
	else {

		theFileStruct = theVectFilesIter->second;

		// make sure not close from a last commit
		if (theFileStruct->theFile == 0) {

			if (vectSize > currentMaxDouble) {
				if (this->resizeDouble(vectSize) < 0) {
					opserr << "FileDatastore::sendVector() - failed in resizeInt()\n";
					return -1;
				}
			}
			char* fileName = new char[strlen(dataBase) + 21];
			static char intName[20];
			strcpy(fileName, dataBase);
			sprintf(intName, "%d.%d", vectSize, commitTag);
			strcat(fileName, ".VECs.");
			strcat(fileName, intName);

			if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
				opserr << "FileDatastore::sendVector() - could not open file\n";
				delete[] fileName;
				return -1;
			}

			delete[] fileName;
		}
	}

	//
	// find location in file to place the data
	//

	fstream* theStream = theFileStruct->theFile;
	bool found = false;
	STREAM_POSITION_TYPE pos = theStream->tellg();
	STREAM_POSITION_TYPE fileEnd = theFileStruct->fileEnd;

	// we first check if the data can go at the end of the file
	// true if commitTag larger than any we have encountered so far
	found = false;


	if (theFileStruct->maxDbTag < dataTag) {
		pos = fileEnd;
		found = true;
		theFileStruct->maxDbTag = dataTag;
	}


	// try current location    
	if (pos < fileEnd && found == false) {
		// #ifdef _WIN32
				// must be a bug in the vc compiler! .. we are here already tellg() above!!
		theStream->seekg(pos, ios::beg);
		// #endif
		theStream->read(data, stepSize);
		if (*(theIntData.dbTag) == dataTag) {
			found = true;
		}
	}

	// we have to search from the beginning of the file  
	if (found == false) {
		*(theIntData.dbTag) = -1;
		pos = sizeof(int);
		theStream->seekg(pos, ios::beg);
		while (pos < fileEnd && found == false) {
			theStream->read(data, stepSize);

			if (*(theIntData.dbTag) == dataTag)
				found = true;
			else
				pos += stepSize;
		}

		if (found == false) {
			pos = fileEnd;
		}
	}

	//
	// we now place the data to be sent into our buffer
	//

	*(theDoubleData.dbTag) = dataTag;
	for (int i = 0; i < vectSize; i++)
		theDoubleData.data[i] = theVector(i);

	//
	// we now write the data
	//

	theStream->seekp(pos, ios::beg); // reset so can go write at the end
	theStream->write(data, stepSize);

	// update the size of file if we have added to eof
	if (theFileStruct->fileEnd <= pos)
		theFileStruct->fileEnd += stepSize;

	return 0;
}

/*
PSEUDOCODE / DETAILED PLAN
- Problem: recvVector was copying values element-by-element using theVector(i) = theDoubleData.data[i].
	Observed symptom: received vector is offset (appears to be the vector sent after two [0,0,0] vectors).
- Hypothesis: element-wise assignment via operator() may interact with Vector internals or be slower/unsafe in this context.
- Fix: copy the raw double data block directly into the Vector storage using GetData() + memcpy when possible.
	Fallback: if GetData() returns nullptr (defensive), revert to element-wise assignment.
- Preserve existing file-search/reading logic; only change the copy step near end of recvVector.
- Keep existing _DEBUG logging intact.

Change: replace the loop
		for (int i = 0; i < vectSize; i++)
				theVector(i) = theDoubleData.data[i];
with a block that uses memcpy from theDoubleData.data into theVector.GetData().
*/

int
FileDatastore::recvVector(int dataTag, int commitTag,
	Vector& theVector,
	ChannelAddress* theAddress)
{
	if (currentCommitTag != commitTag)
		this->resetFilePointers();

	currentCommitTag = commitTag;

	FileDatastoreOutputFile* theFileStruct;

	//
	// next we see if we already have this file; 
	//  if not we need to create data structure & open it
	//  if we have data structure, need to check file is opened (we close in a commit)
	//

	// we first check Vector not too big
	int vectSize = theVector.Size();
	int stepSize = sizeof(int) + vectSize * sizeof(double);

	theVectFilesIter = theVectFiles.find(vectSize);
	if (theVectFilesIter == theVectFiles.end()) {

		// we first check if we need to resize recv buffer
		if (vectSize > currentMaxDouble) {
			if (this->resizeDouble(vectSize) < 0) {
				opserr << "FileDatastore::recvVectrix() - failed in resizeDouble()\n";
				return -1;
			}
		}

		char* fileName = new char[strlen(dataBase) + 21];
		theFileStruct = new FileDatastoreOutputFile;

		if (fileName == 0 || theFileStruct == 0) {
			opserr << "FileDatastore::recvVectrix() - out of memory\n";
			return -1;
		}

		static char intName[20];
		strcpy(fileName, dataBase);
		sprintf(intName, "%d.%d", vectSize, commitTag);
		strcat(fileName, ".VECs.");
		strcat(fileName, intName);

		if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
			opserr << "FileDatastore::recvVectrix() - could not open file\n";
			delete[] fileName;
			return -1;
		}
		else
			theVectFiles.insert(MAP_FILES_TYPE(vectSize, theFileStruct));

		delete[] fileName;

	}
	else {

		theFileStruct = theVectFilesIter->second;

		// make sure not close from a last commit
		if (theFileStruct->theFile == 0) {

			if (vectSize > currentMaxDouble) {
				if (this->resizeDouble(vectSize) < 0) {
					opserr << "FileDatastore::recvVectrix() - failed in resizeInt()\n";
					return -1;
				}
			}
			char* fileName = new char[strlen(dataBase) + 21];
			static char intName[20];
			strcpy(fileName, dataBase);
			sprintf(intName, "%d.%d", vectSize, commitTag);
			strcat(fileName, ".VECs.");
			strcat(fileName, intName);

			if (this->openFile(fileName, theFileStruct, stepSize) < 0) {
				opserr << "FileDatastore::recvVectrix() - could not open file\n";
				delete[] fileName;
				return -1;
			}

			delete[] fileName;
		}
	}

	fstream* theStream = theFileStruct->theFile;
	STREAM_POSITION_TYPE fileEnd = theFileStruct->fileEnd;
	STREAM_POSITION_TYPE pos = theStream->tellg();

	//
	// find location in file to place the data
	//

	bool found = false;

	// we try the current file position first

	if (pos < fileEnd) {
		theStream->read(data, stepSize);
		if ((*(theIntData.dbTag) == dataTag)) {
			found = true;
			pos += stepSize;
		}
	}

	// we must search from the beginning of the file
	if (found == false) {
		pos = sizeof(int);
		theStream->seekg(pos, ios::beg);
		while ((pos < fileEnd) && (found == false)) {
			theStream->read(data, stepSize);
			if (*(theIntData.dbTag) == dataTag)
				found = true;
			else
				pos += stepSize;
		}
	}

	//
	// find location in file to place the data
	//

	/*bool found = false;
	STREAM_POSITION_TYPE pos = theStream->tellg();

	// Fast path: try the current get position first (if it's within range).
	if (pos < fileEnd) {
		theStream->seekg(pos, ios::beg);
		theStream->read(data, stepSize);
		if (theStream->good() && (*(theIntData.dbTag) == dataTag)) {
			found = true;
		}
	}

	// Robust fallback: explicitly seek to each record position and read,
	// so we don't rely on the stream's internal cursor state.
	if (!found) {
		for (STREAM_POSITION_TYPE p = sizeof(int); p + stepSize <= fileEnd; p += stepSize) {
			theStream->seekg(p, ios::beg);
			theStream->read(data, stepSize);
			if (!theStream->good()) {
				theStream->clear(); // clear flags and continue (defensive)
				continue;
			}
			if (*(theIntData.dbTag) == dataTag) {
				found = true;
				pos = p;
				break;
			}
		}
	}*/
	if (found == false) {
		opserr << "FileDatastore::recvVector() - failed\n";
		return -1;
	}

	// Copy the retrieved raw doubles into theVector efficiently and safely.
	// Prefer a block copy into the Vector's internal storage; fallback to element-wise assign.
	double* vecData = theVector.GetData();
	if (vecData != 0) {
		// memcpy the raw doubles; theDoubleData.data points into the internal data buffer
		// (ensure we copy at most vectSize elements)
		std::memcpy(vecData, theDoubleData.data, static_cast<size_t>(vectSize) * sizeof(double));
	}
	else {
		// defensive fallback
		for (int i = 0; i < vectSize; i++)
			theVector(i) = theDoubleData.data[i];
	}

#ifdef _DEBUG
	std::ostringstream oss;
	oss << "Vector: ";
	writeVectorToStream(oss, theVector);
	appendLogEntry(dataTag, dataBase, oss.str());
#endif
	return 0;
}


using namespace std;
int
FileDatastore::createTable(const char* tableName, int numColumns, char* columns[])
{
	// open the file
	int res = 0;
	char* fileName = new char[strlen(tableName) + strlen(dataBase) + 10];
	if (fileName == 0) {
		opserr << "FileDatastore::insertData - out of memory; failed to open file: " << fileName << endln;
		return -1;
	}

	strcpy(fileName, dataBase);
	strcat(fileName, ".");
	strcat(fileName, tableName);

	ofstream table;
	table.open(fileName, ios::out | ios::trunc);

	if (table.bad() == true || table.is_open() == false) {
		opserr << "FileDatastore::insertData - failed to open file: " << fileName << endln;
		delete[] fileName;
		res = -1;

	}

	// write the data
	for (int i = 0; i < numColumns; i++) {
		table << columns[i] << "\t";
	}
	table << "\n";
	table.close();

	delete[] fileName;

	return res;
}

int
FileDatastore::insertData(const char* tableName, char* columns[],
	int commitTag, const Vector& data)
{
	// open the file
	char* fileName = new char[strlen(tableName) + strlen(dataBase) + 10];
	if (fileName == 0) {
		opserr << "FileDatastore::insertData - out of memory; failed to open file: " << fileName << endln;
		return -1;
	}

	strcpy(fileName, dataBase);
	strcat(fileName, ".");
	strcat(fileName, tableName);

	ofstream table;
	table.open(fileName, ios::app);
	if (table.bad() == true || table.is_open() == false) {
		opserr << "FileDatastore::insertData - failed to open file: " << fileName << endln;
		delete[] fileName;
		return -1;
	}

	table << setiosflags(ios::scientific);
	table << std::setprecision(16);

	// write the data
	for (int i = 0; i < data.Size(); i++) {
		table << data(i) << "\t";
	}

	table << "\n";
	table.close();

	delete[] fileName;
	return 0;
}


int
FileDatastore::getData(const char* tableName, char* columns[], int commitTag, Vector& data)
{
	return 0;
}



/*******************************************************************
 *              MISC METHODS & FUNCTIONS FOR OPENING THE FILE      *
 *******************************************************************/

int
FileDatastore::openFile(char* fileName, FileDatastoreOutputFile* theFileStruct, int dataSize)
{
	fstream* res = new fstream();
	if (res == 0) {
		opserr << "FileDatastore::openFile - out of memory; failed to open file: " << fileName << endln;
		return 0;
	}

	res->open(fileName, ios::in | ios::out | ios::binary);

	// if file did not exist, need to pass trunc flag to open it
	if (res->bad() == true || res->is_open() == false) {
		// delete & new again for unix gcc compiler to work!
		delete res;
		res = new fstream();
		if (res == 0) {
			opserr << "FileDatastore::openFile - out of memory; failed to open file: " << fileName << endln;
			theFileStruct->theFile = res;
			return -1;
		}
		res->open(fileName, ios::in | ios::out | ios::trunc | ios::binary);
	}

	if (res->bad() == true || res->is_open() == false) {
		opserr << "FATAL - FileDatastore::openFile() - could not open file " << fileName << endln;
		delete res;
		theFileStruct->theFile = 0;
		return -1;
	}

	// set the position for writing to eof
	res->seekp(0, ios::end);
	STREAM_POSITION_TYPE fileEnd = res->tellp();
	int maxDataTag = 0;

	if (fileEnd == 0 || fileEnd == -1) {
		*(theIntData.dbTag) = maxDataTag;
		res->write(data, sizeof(int));
		fileEnd = sizeof(int);
		maxDataTag = -1;
	}
	else {
		res->seekg(0, ios::beg);
		res->read(data, sizeof(int));
		maxDataTag = *(theIntData.dbTag);
	}

	// move to start of data part
	res->seekp(sizeof(int), ios::beg);
	res->seekg(sizeof(int), ios::beg);

	// fill in the structure data
	theFileStruct->theFile = res;
	theFileStruct->fileEnd = fileEnd;

	theFileStruct->maxDbTag = maxDataTag;

	return 0;
}

int
FileDatastore::resizeInt(int newSize) {
	int sizeOfChar = sizeof(char);
	int sizeOfInt = sizeof(int);
	int sizeOfDouble = sizeof(double);
	newSize = (newSize + 1) * sizeOfInt / sizeOfChar;

	if (newSize < sizeData)
		return 0; // already big enough

	if (newSize <= 0) {
		opserr << "FileDatastore::resizeInt(int newSize) - invalidSize " << newSize << endln;
		return -1; // invalid size
	}

	if (data != 0)
		delete[] data;

	data = new char[newSize];
	if (data == 0) {
		opserr << "FileDatastore::resizeInt(int newSize) - out of memory for size: " << newSize << endln;
		return -1;
	}

	sizeData = newSize;

	currentMaxInt = (sizeData / sizeOfChar - sizeOfInt) / sizeOfInt;
	currentMaxDouble = (sizeData / sizeOfChar - sizeOfInt) / sizeOfDouble;
	char* dataPtr = &data[sizeof(int)];
	theIntData.dbTag = (int*)data;
	theIntData.data = (int*)dataPtr;
	theDoubleData.dbTag = (int*)data;
	theDoubleData.data = (double*)dataPtr;

	return 0;
}

int
FileDatastore::resizeDouble(int newSize) {
	int sizeOfChar = sizeof(char);
	int sizeOfInt = sizeof(int);
	int sizeOfDouble = sizeof(double);
	newSize = (newSize * sizeOfDouble + sizeOfInt) / sizeOfChar;

	if (newSize < sizeData)
		return 0; // already big enough

	if (newSize <= 0) {
		opserr << "FileDatastore::resizeInt(int newSize) - invalidSize " << newSize << endln;
		return -1; // invalid size
	}

	if (data != 0)
		delete[] data;

	data = new char[newSize];
	if (data == 0) {
		opserr << "FileDatastore::resizeInt(int newSize) - out of memory for size: " << newSize << endln;
		return -1;
	}

	sizeData = newSize;

	currentMaxInt = (sizeOfChar * sizeData - sizeOfInt) / sizeOfInt;
	currentMaxDouble = (sizeOfChar * sizeData - sizeOfInt) / sizeOfDouble;
	char* dataPtr = &data[sizeof(int)];
	theIntData.dbTag = (int*)data;
	theIntData.data = (int*)(dataPtr);
	theDoubleData.dbTag = (int*)data;
	theDoubleData.data = (double*)(dataPtr);
	return 0;
}

#undef _DEBUG