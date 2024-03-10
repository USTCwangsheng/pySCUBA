/*
 * blkselectionparser.h
 *
 *  Created on: 2019年12月12日
 *      Author: hyiu
 */

#ifndef IBLOCK_BLKSELECTIONPARSER_H_
#define IBLOCK_BLKSELECTIONPARSER_H_
#include <string>
#include <map>
#include <set>
#include <vector>
namespace NSPintrct{
class IntrctMol;
/**
 * defines a selected set of blocks in the molecular system
 */
struct BlkSelector{
	typedef std::map<int,std::set<int>> SelectedBlks; //<* the key is the chain  id, the value is the selected residue in that chain
	BlkSelector(){;}
	BlkSelector(const std::string &selectionstr){
		addselection(selectionstr);
	}
	SelectedBlks selectedblks; //<* selected blocks
	void clear(){selectedblks.clear();}
/**
 * parse a string to obtain selected positions in chains.
 * for example, for an input string "chain0 5, 7-10, chain2 3,5,7,9-10",
 * the resulting map will contain:
 * {{0,{5,7,8,9,10}},{2,{3,5,7,9,10}}}*/
	void addselection(const std::string & selectionstr);
	std::vector<int> selectatoms(const IntrctMol &imol) const;
	std::vector<int> selectmcatoms(const IntrctMol &imol) const;
	std::vector<int> selectscatoms(const IntrctMol &imol) const;
};
}



#endif /* IBLOCK_BLKSELECTIONPARSER_H_ */
