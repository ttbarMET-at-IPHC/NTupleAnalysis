#include "EventReco/interface/TablesMaker.h"
#include <cstring>

TablesMaker::TablesMaker() {}
TablesMaker::~TablesMaker() {}

void TablesMaker::Add(const std::vector<std::string> headLine)
{
	tableEntries = headLine;
}

void TablesMaker::Add(const std::string chars,const std::vector<double> entry)
{
	lineNames.push_back(chars);
	entries.push_back(entry);
}

void TablesMaker::PrintTable()
{
	tableSize_ = tableEntries.size()+1.;
	for (int i=0; i<tableSize_; i++)
	{
		if (i==0)
		{
			for (int j=0; j<20. ; j++) std::cout << " ";
		}
		else if (i== (tableSize_-1.)) 
		{
			std::cout << tableEntries[i-1] << std::endl; 
		}
		else
		{
			std::cout <<tableEntries[i-1];
			for (int j=0; j<(20.-tableEntries[i-1].size()) ; j++) std::cout << " ";
		}
	}
	nOflines = lineNames.size();
	for (int j=0;j<nOflines;j++)
	{
		for (int i=0;i<tableSize_;i++)
		{
			if (i==0)
			{
				std::cout << lineNames[j];
				nbOfSpaces = 20. - lineNames[j].size();
				for (int k=0;k<nbOfSpaces;k++) std::cout <<" ";
			}
			else if (i==tableSize_-1)
			{
				for(int k=0;k<(20.-sizeof((entries[j])[i-1]));k++) std::cout <<" ";	
				std::cout << (entries[j])[i-1]<< std::endl;
			}
			else
			{	
				std::cout <<(entries[j])[i-1];
			}
		}
	}
}
