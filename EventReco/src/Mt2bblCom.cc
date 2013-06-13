#include "EventReco/interface/Mt2bblCom.h"


Mt2bblCom::Mt2bblCom(const std::vector<IPHCTree::NTJet>& VBJet,
										 const std::vector<IPHCTree::NTJet>& VJet,
									 	 const std::vector<IPHCTree::NTElectron>& VElectron,
										 const std::vector<IPHCTree::NTMuon>& VMuon)
{
	NTBJet=VBJet;
	NTJet=VJet;
	NTElectron=VElectron;
	NTMuon=VMuon;
}

Mt2bblCom::~Mt2bblCom() {}


//Calcul des Mt2 bl et b
void Mt2bblCom::ComputeMt2bbl()
{
					Mt2bCalc.clear();
					Mt2blCalc.clear();
					ctr=0;
					it=0;
					// UN SEUL B-JET
					if (NTBJet.size()==1)
					{
						if (NTElectron.size()==0)
						{
							while(ctr<=2)
							{
								if (NTBJet[0].p4!=NTJet[it].p4)
								{
									Mt2Com::SetVisA(NTBJet[0].p4 + NTMuon[0].p4);
									Mt2Com::SetVisB(NTJet[it].p4);
									Mt2Com::SetMInvisMass(0);
									Mt2blCalc.push_back(Mt2Com::GetMt2_332());
	
									Mt2Com::SetVisA(NTJet[it].p4 + NTMuon[0].p4);
									Mt2Com::SetVisB(NTBJet[0].p4);
									Mt2Com::SetMInvisMass(0);
									Mt2blCalc.push_back(Mt2Com::GetMt2_332());

									Mt2Com::SetVisA(NTBJet[0].p4);
									Mt2Com::SetVisB(NTJet[it].p4);
									Mt2Com::SetMInvisMass(0);
									Mt2bCalc.push_back(Mt2Com::GetMt2_332());
				
									it ++;
									ctr ++;
								}
								else it ++;
							}

							//	cout << "Valeur de Mt2 1 B-jet= " << *std::min_element(Mt2Calc.begin(),Mt2Calc.end())<< endl;
						}
						else if (NTMuon.size()==0)
						{
							while(ctr<=2)
							{
								if (NTBJet[0].p4!=NTJet[it].p4)
								{
									Mt2Com::SetVisA(NTBJet[0].p4 + NTElectron[0].p4);
									Mt2Com::SetVisB(NTJet[it].p4);
									Mt2Com::SetMInvisMass(0);
									Mt2blCalc.push_back(Mt2Com::GetMt2_332());
	
									Mt2Com::SetVisA(NTJet[it].p4 + NTElectron[0].p4);
									Mt2Com::SetVisB(NTBJet[0].p4);
									Mt2Com::SetMInvisMass(0);
									Mt2blCalc.push_back(Mt2Com::GetMt2_332());

									Mt2Com::SetVisA(NTBJet[0].p4);
									Mt2Com::SetVisB(NTJet[it].p4);
									Mt2Com::SetMInvisMass(0);
									Mt2bCalc.push_back(Mt2Com::GetMt2_332());
				
									it ++;
									ctr ++;
								}
								else it ++;
							}
							
							//	cout << "Valeur de Mt2 1 B-jet= " << .min_element(Mt2Calc.begin(),Mt2Calc.end())<< endl;
						}
						else cout << "Erreur : pas de lepton dans l'evenement" << endl;
					}
					// DEUX B-JETS					
					else if (NTBJet.size()==2)
					{
						if (NTElectron.size()==0)
						{
							Mt2Com::SetVisA(NTBJet[0].p4 + NTMuon[0].p4);
							Mt2Com::SetVisB(NTBJet[1].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());
	
							Mt2Com::SetVisA(NTBJet[1].p4 + NTMuon[0].p4);
							Mt2Com::SetVisB(NTBJet[0].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTBJet[0].p4);
							Mt2Com::SetVisB(NTBJet[1].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2bCalc.push_back(Mt2Com::GetMt2_332());

							//	cout << "Valeur de Mt2 2 B-jets= " << *std::min_element(Mt2Calc.begin(),Mt2Calc.end())<< endl;
						}
						else if (NTMuon.size()==0)
						{
							Mt2Com::SetVisA(NTBJet[0].p4 + NTElectron[0].p4);
							Mt2Com::SetVisB(NTBJet[1].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());
	
							Mt2Com::SetVisA(NTBJet[1].p4 + NTElectron[0].p4);
							Mt2Com::SetVisB(NTBJet[0].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTBJet[0].p4);
							Mt2Com::SetVisB(NTBJet[1].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2bCalc.push_back(Mt2Com::GetMt2_332());

							//	cout << "Valeur de Mt2 2 B-jets = " << *std::min_element(Mt2Calc.begin(),Mt2Calc.end())<< endl;
						}
						else cout << "Erreur : pas de lepton dans l'evenement" << endl;
					}
					// TROIS B-JETS OU PLUS				
					else
					{
						if (NTElectron.size()==0)
						{
							Mt2Com::SetVisA(NTJet[0].p4 + NTMuon[0].p4);
							Mt2Com::SetVisB(NTJet[1].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTJet[1].p4 + NTMuon[0].p4);
							Mt2Com::SetVisB(NTJet[0].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTJet[0].p4 + NTMuon[0].p4);
							Mt2Com::SetVisB(NTJet[2].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTJet[2].p4 + NTMuon[0].p4);
							Mt2Com::SetVisB(NTJet[0].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTJet[1].p4 + NTMuon[0].p4);
							Mt2Com::SetVisB(NTJet[2].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTJet[2].p4 + NTMuon[0].p4);
							Mt2Com::SetVisB(NTJet[1].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTBJet[0].p4);
							Mt2Com::SetVisB(NTBJet[1].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2bCalc.push_back(Mt2Com::GetMt2_332());
								
							Mt2Com::SetVisA(NTBJet[1].p4);
							Mt2Com::SetVisB(NTBJet[2].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2bCalc.push_back(Mt2Com::GetMt2_332());
							
							Mt2Com::SetVisA(NTBJet[0].p4);
							Mt2Com::SetVisB(NTBJet[2].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2bCalc.push_back(Mt2Com::GetMt2_332());

							//	cout << "Valeur de Mt2Com::3 B-jets ou plus= " << *std::min_element(Mt2Calc.begin(),Mt2Calc.end())<< endl;
						}
						else if (NTMuon.size()==0)
						{
							Mt2Com::SetVisA(NTJet[0].p4 + NTElectron[0].p4);
							Mt2Com::SetVisB(NTJet[1].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTJet[1].p4 + NTElectron[0].p4);
							Mt2Com::SetVisB(NTJet[0].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTJet[0].p4 + NTElectron[0].p4);
							Mt2Com::SetVisB(NTJet[2].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTJet[2].p4 + NTElectron[0].p4);
							Mt2Com::SetVisB(NTJet[0].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTJet[1].p4 + NTElectron[0].p4);
							Mt2Com::SetVisB(NTJet[2].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());
	
							Mt2Com::SetVisA(NTJet[2].p4 + NTElectron[0].p4);
							Mt2Com::SetVisB(NTJet[1].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2blCalc.push_back(Mt2Com::GetMt2_332());

							Mt2Com::SetVisA(NTBJet[0].p4);
							Mt2Com::SetVisB(NTBJet[1].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2bCalc.push_back(Mt2Com::GetMt2_332());
								
							Mt2Com::SetVisA(NTBJet[1].p4);
							Mt2Com::SetVisB(NTBJet[2].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2bCalc.push_back(Mt2Com::GetMt2_332());
							
							Mt2Com::SetVisA(NTBJet[0].p4);
							Mt2Com::SetVisB(NTBJet[2].p4);
							Mt2Com::SetMInvisMass(0);
							Mt2bCalc.push_back(Mt2Com::GetMt2_332());

						//		cout << "Valeur de Mt2 3 B-jets ou plus= " << *std::min_element(Mt2Calc.begin(),Mt2Calc.end())<< endl;
						} 
						else std::cout << "Erreur : pas de lepton dans l'evenement" << std::endl;
					}
}

double Mt2bblCom::GetMt2bl()
{
	Mt2bblCom::ComputeMt2bbl();
	return *std::min_element(Mt2blCalc.begin(), Mt2blCalc.end());
}

double Mt2bblCom::GetMt2b()
{
	Mt2bblCom::ComputeMt2bbl();
	return *std::min_element(Mt2bCalc.begin(), Mt2bCalc.end());
}
