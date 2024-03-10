/*
 * molpara.h
 *
 *  Created on: 2019年12月17日
 *      Author: hyiu
 */

#ifndef IBLOCK_MOLSYSTMPARA_H_
#define IBLOCK_MOLSYSTMPARA_H_
#include "dstl/vec2d.h"
#include "iblock/blkselector.h"
#include "dataio/controlfile.h"
#include <vector>
#include <set>
#include <string>
namespace NSPintrct{
/**
 * Parameters defining the set up of the molecular system
 *
 * The fixed parts of the molecular system can be assigned by using the
 * fixedmainchain residues and fixedresidues selections. The remaining
 * unselected part will be active.
 *
 * The active parts can be assigned by using the activeresidues and
 * sidechainactive residues selections. The remaining unselected part will
 * be fixed.
 *
 * The fixed parts and active parts CANNOT be specified simultaneously.
 */
struct MolSystmPara{
	std::vector<std::string> sequence1l;
	std::string pdbstart{""};  //*< the starting PDB file. Seqeunces and coordinates may come from this file.
	bool changepdbsequence{false};
	std::string jobname{""};
//	BlkSelector::SelectedBlks cispositions;
	BlkSelector::SelectedBlks fixedmainchainresidues; //*<main chain of some resdiues may be fixed
	BlkSelector::SelectedBlks fixedresidues; //*<some residues may be completely fixed
	BlkSelector::SelectedBlks activeresidues; //*< all atoms in some residues can move
	BlkSelector::SelectedBlks sidechainactiveresidues; //*< sidechain atoms in some residues can move
	BlkSelector::SelectedBlks softsidechainresidues; //*< residues use soft sidechain packing interactions
	MolSystmPara(){;}
	MolSystmPara(const std::vector<std::string> & controllines);
};

/**
 * read the parameters from an input control file
 */
inline MolSystmPara makemolsystmparam(const std::string &jobname,const NSPdataio::ControlFile &cf,
		const std::string & controlname="MolSystmPar"){
		std::vector<std::string> lines=cf.getcontrolines(controlname);
		if (lines.empty()) {
			std::cout << "NO parameters in MolSystmPar" << std::endl;
			exit(1);
		}
		if(!jobname.empty()) lines.push_back("JobName = "+jobname);
		return MolSystmPara(lines);
}
}

#endif /* IBLOCK_MOLSYSTMPARA_H_ */
