/*
 * pdbatomrecord.h
 *
 *  Created on: 2016年11月16日
 *      Author: hyliu
 */

#ifndef PROTEINREP_PDBRECORD_H_
#define PROTEINREP_PDBRECORD_H_
#include <string>
#include <vector>
namespace NSPproteinrep{
#ifndef NOID
#define NOID -100000
#endif
#ifndef CHAINID_DEFAULT
#define CHAINID_DEFAULT 'A'
#endif
/*! data in an ATOM or HETATM record
 *
 */
struct PdbRecord {
	/*!fields in one line
	 *
	 */
	enum Field {
		LABEL,        /*!< enum value LABEL, whose position is column 1-6 in pdb file.*/
		ATOMID,       /*!< enum value ATOMID, whose position is column 7-12 in pdb file.*/
		NAMESYMBOL,   /*!< enum value NAMESYMBOL, whose position is column 13-14 in pdb file.*/
		NAMEMODIFIER, /*!< enum value NAMEMODIFIER, whose position is column 15-16 in pdb file.*/
		CONFORMERID,  /*!< enum value CONFORMERID, whose position is column 17 in pdb file.*/
		RESIDUENAME,  /*!< enum value RESIDUENAME, whose position is column 18-21 in pdb file.*/
		CHAINID,      /*!< enum value CHAINID, whose position is column 22 in pdb file.*/
		RESIDUEID,    /*!< enum value RESIDUEID, whose position is column 23-26 in pdb file.*/
		INSERTIONID,  /*!< enum value INSERTIONID, whose position is column 27 in pdb file.*/
		X,            /*!< enum value X, whose position is column 28-38 in pdb file.*/
		Y,            /*!< enum value Y, whose position is column 39-46 in pdb file.*/
		Z,            /*!< enum value Z, whose position is column 47-54 in pdb file.*/
		OCCUPATION,   /*!< enum value OCCUPATION, whose position is column 55-60 in pdb file.*/
		BFACTOR,      /*!< enum value BFACTOR, temperature factor,whose position is column 61-66 in pdb file.*/
		SEGMENT,      /*!< enum value SEGMENT, whose position is column 67-76 in pdb file.*/
		ELEMENTNAME   /*!< enum value ELEMENTNAME, element symbol, whose position is column 77-78 in pdb file.*/
	};
	/*!default constructor
	 *
	 */
	PdbRecord() {
		;
	}
	PdbRecord(const std::string & line);
    /*! initialize a PdbRecord instance
     *
     * @param fields a vector contains fields
     */
	void init(const std::vector<std::string> &fields);
	std::string toString() const;
    // member variables
	std::string label { "" };
	std::string atomname { "" };  //atomname=namesymbol+namemodifier
	std::string namesymbol { "" };
	std::string namemodifier { "" };
	std::string residuename { "" };
	int atomid { NOID };
	char chainid { ' ' };
	int residueid { NOID };
	char conformerid { ' ' };
	char insertionid { ' ' };
	double x { 0.0 }, y { 0.0 }, z { 0.0 };
	double occupation { 1.0 };
	double bfactor { 0.0 };
	char elementname[2] { ' ', 'X' };
};
/*! make a pdb record with given atom information
 *
 * @tparam ATOMKEY
 * @tparam XYZ
 * @param key      contains position number, atomname, chain number(default 0) , residue type(default "ANY"),
 *                 Rotamer Code Number(default 0)
 * @param xyz      vector stored coordinates of xyz
 * @param atomid   atom id
 * @param residueidshift
 * @param elementsymbolsize  size of element symbol
 * @return  character string identifying the pdb recode
 */
template<typename ATOMKEY,typename XYZ>
PdbRecord make_pdbrecord(const typename ATOMKEY::Key & key, const XYZ & xyz,int atomid=0, int residueidshift=0,
		int elementsymbolsize=1) {
	PdbRecord record;
	const std::vector<char> chainids{'A','B','C','D','E','F','G','H','I','J'};
	record.label="ATOM";
	record.atomname=ATOMKEY::atomName(key);     //get atomname from key, and assign the atomname to record.atomname
	record.namesymbol=record.atomname.substr(0,elementsymbolsize);
	if(elementsymbolsize > 1) {
		record.elementname[0]= record.namesymbol[0];
		record.elementname[1]=record.namesymbol[1];
	} else {
		record.elementname[1]=record.namesymbol[0];
	}
	record.namemodifier=record.atomname.substr(elementsymbolsize);
	record.residuename=ATOMKEY::residueName(key);
	record.atomid=atomid;
	record.residueid = ATOMKEY::posiNumber(key) + residueidshift;
	record.x=xyz[0];
	record.y=xyz[1];
	record.z=xyz[2];
	int chainid=ATOMKEY::chainNumber(key);
	if(chainid <10) record.chainid=chainids[chainid];
	return PdbRecord(record.toString());
}
}

#endif /* PROTEINREP_PDBRECORD_H_ */
