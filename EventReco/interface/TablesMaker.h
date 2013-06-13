#ifndef TablesMaker_H
#define TablesMaker_H

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <vector>

class TablesMaker
{
	public :

	//Constructor
	TablesMaker();

	//Destructor
	~TablesMaker();

	//Method members

	void Add(const std::vector<std::string> headLine);
	void Add(const std::string chars, const std::vector<double> entry);
	void PrintTable ();

	private :
	int nbOfSpaces;
	std::vector<std::string> tableEntries;
	int tableSize_;
	int nOflines;
	std::vector< std::vector<double> > entries;
	std::vector<std::string> lineNames;	
};
#endif
