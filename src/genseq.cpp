/*
    loresim
    Copyright (C) 2015 German Tischler

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/random/GaussianRandom.hpp>
#include <libmaus2/fastx/acgtnMap.hpp>
#include <set>

enum err_type_enum
{
    err_subst,
    err_ins,
    err_del
};

static bool contains(std::vector<err_type_enum> const & V, err_type_enum type)
{
    for ( uint64_t i = 0; i < V.size(); ++i )
        if ( V[i] == type )
            return true;
    return false;
}

enum state_enum
{
    state_error_low = 0,
    state_error_high = 1
};


int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		
		uint64_t const blocksizeavg = arginfo.getValueUnsignedNumeric<uint64_t>("blocksizeavg",8*1024);
		uint64_t const blocksizestd = arginfo.getValueUnsignedNumeric<uint64_t>("blocksizestd",4*1024);
		uint64_t const blockmin = arginfo.getValueUnsignedNumeric<uint64_t>("blockmin",512);
		double repfrag = arginfo.getValue<double>("repfrag",0.01);	
		double reprepfrag = arginfo.getValue<double>("reprepfrag",0.5);
		uint64_t const numblocks = arginfo.getValueUnsignedNumeric<uint64_t>("numblocks",1024);
		double const eratesq = arginfo.getValue<double>("eratesq",0.001);
		double const erateavg = arginfo.getValue<double>("erateavg",0.001);
		double substrate = arginfo.getValue<double>("substrate",.2);
		double delrate = arginfo.getValue<double>("delrate",.3);
		double insrate = arginfo.getValue<double>("insrate",.5);
		
		std::vector<uint64_t> blockids;
		std::vector<std::string> blocks;
		std::set<uint64_t> repblockset;
		libmaus2::random::Random::setup();
		
		for ( uint64_t i = 0; i < numblocks; ++i )
		{
			if ( blocks.size() && libmaus2::random::UniformUnitRandom::uniformUnitRandom() < repfrag )
			{
				uint64_t repblock;
				
				if ( repblockset.size() && libmaus2::random::UniformUnitRandom::uniformUnitRandom() < reprepfrag )
				{
					uint64_t repoff = libmaus2::random::Random::rand64() % repblockset.size();
					std::set<uint64_t>::const_iterator ita = repblockset.begin();
					while ( repoff-- )
						++ita;
					repblock = *ita;
				}
				else
				{
					repblock = libmaus2::random::Random::rand64() % blocks.size();
				}

				std::cerr << "[D] rep " << repblock << " of size " << blocks[repblock].size() << std::endl;
				blockids.push_back(repblock);				
				repblockset.insert(repblock);
			}
			else
			{
				double dlength = -1;

				while ( dlength <= blockmin )
					dlength = libmaus2::random::GaussianRandom::random(blocksizestd,blocksizeavg);

				uint64_t const length = dlength;
				
				std::vector<char> V(length);
				for ( uint64_t j = 0; j < V.size(); ++j )
					V[j] = libmaus2::fastx::remapChar(libmaus2::random::Random::rand8()%4);
					
				std::cerr << "[D] new block of size " << V.size() << std::endl;
				
				blockids.push_back(blocks.size());
				blocks.push_back(std::string(V.begin(),V.end()));
			}
		}
		
		std::cerr << "[D] all blocks generated" << std::endl;
		
		for ( uint64_t i = 0; i < blockids.size(); ++i )
		{
			
			std::string sub = blocks[blockids[i]];
			
			std::map<uint64_t,std::vector<err_type_enum> > errM;
			for ( uint64_t i = 0; i < sub.size(); ++i )
			    errM[i] = std::vector<err_type_enum>(0);

			std::map<int,double> err_prob_cumul;
			double ratesum = substrate + delrate + insrate;
			substrate /= ratesum;
			delrate /= ratesum;
			insrate /= ratesum;
        
			err_prob_cumul[err_subst] = substrate;
			err_prob_cumul[err_del] = substrate + delrate;
			err_prob_cumul[err_ins] = substrate + delrate + insrate;

			double erate = -1;
			while ( (erate = libmaus2::random::GaussianRandom::random(eratesq,erateavg)) < 0 )
			{    
			}
			    
			uint64_t const numerr = std::floor(erate * sub.size() + 0.5);

			uint64_t errplaced = 0;
			    
			while ( errplaced < numerr )
			{
				double const p = libmaus2::random::UniformUnitRandom::uniformUnitRandom();
				uint64_t const errpos = std::floor(libmaus2::random::UniformUnitRandom::uniformUnitRandom() * (sub.size()-1) + 0.5);
				
				if ( p < err_prob_cumul[err_subst]  )
				{
					// no error yet?
					if ( errM.find(errpos) == errM.end() )
					{
						errM[errpos].push_back(err_subst);
						errplaced++;
					}
					// no error so far a deletion or substitution?
					else if ( 
						!contains(errM.find(errpos)->second,err_del)
						&&
						!contains(errM.find(errpos)->second,err_subst)
					)
					{
						errM[errpos].push_back(err_subst);
						errplaced++;                        
					}
				}
				else if ( p < err_prob_cumul[err_del] )
				{
					// no error yet?
					if ( errM.find(errpos) == errM.end() )
					{
						errM[errpos].push_back(err_del);
						errplaced++;
					}                
					// no error so far a deletion or substitution?
					else if ( 
						!contains(errM.find(errpos)->second,err_del)
						&&
						!contains(errM.find(errpos)->second,err_subst)
					)
					{
						errM[errpos].push_back(err_subst);
						errplaced++;                        
					}
				}
				else
				{
					errM[errpos].push_back(err_ins);
					errplaced++;
				}
			}

			std::ostringstream baseostr;
			for ( std::map<uint64_t,std::vector<err_type_enum> >::const_iterator ita = errM.begin(); ita != errM.end(); ++ita )
			{
				uint64_t const pos = ita->first;
				std::vector<err_type_enum> const & errV = ita->second;
                    
				// insertions                
				for ( uint64_t i = 0; i < errV.size(); ++i )
					if ( errV[i] == err_ins )
					{
						int const v = libmaus2::random::Random::rand8() & 3;
						int insbase = -1;
						switch ( v )
						{
							case 0: insbase = 'A'; break;
							case 1: insbase = 'C'; break;
							case 2: insbase = 'G'; break;
							case 3: insbase = 'T'; break;
						}
						baseostr.put(insbase);
					}

					// base not deleted?
					if ( !contains(errV,err_del) )
					{
						if ( contains(errV,err_subst) )
						{
							int const v = libmaus2::random::Random::rand8() & 3;
							int insbase = -1;
							switch ( v )
							{
								case 0: insbase = 'A'; break;
								case 1: insbase = 'C'; break;
								case 2: insbase = 'G'; break;
								case 3: insbase = 'T'; break;
							}                    
							baseostr.put(insbase);
						}
						else
						{
							baseostr.put(sub[pos]); 
						}
					}
					else
					{
					}
			}
			    
			std::cout << ">seq_" << blockids[i] << "_" << baseostr.str().size() << "\n";
			std::cout << baseostr.str() << "\n";	
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		throw ex;
	}
}
