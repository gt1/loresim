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
#include <libmaus2/fastx/StreamFastAReader.hpp>

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

int main(int argc, char ** argv)
{
    try
    {
        libmaus2::util::ArgInfo const arginfo(argc,argv);
        uint64_t randseed = arginfo.getValueUnsignedNumeric<uint64_t>("randomseed",time(0));
        libmaus2::random::Random::setup(randseed);

        std::map<int,double> err_prob_cumul;
        
        double substrate = arginfo.getValue<double>("substrate",.2);
        double delrate = arginfo.getValue<double>("delrate",.3);
        double insrate = arginfo.getValue<double>("insrate",.5);
        double ratesum = substrate + delrate + insrate;
        substrate /= ratesum;
        delrate /= ratesum;
        insrate /= ratesum;
        
        err_prob_cumul[err_subst] = substrate;
        err_prob_cumul[err_del] = substrate + delrate;
        err_prob_cumul[err_ins] = substrate + delrate + insrate;
        
        if ( std::abs ( err_prob_cumul.find(err_ins)->second - 1.0 ) > 1e-3 )
        {
            libmaus2::exception::LibMausException lme;
            lme.getStream() << "subst, del and insrate cannot be normalised to sum 1" << std::endl;
            lme.finish();
            throw lme;
        }

        double const droprate = arginfo.getValue<double>("droprate",0.01);        
        uint64_t const numtraversals = arginfo.getValueUnsignedNumeric<uint64_t>("numtraversals",1);
        double const readlenavg = arginfo.getValue<double>("readlenavg", 15000);
        double const readlenstddev = arginfo.getValue<double>("readlenstddev", 3000);
        double const keeplowstate = arginfo.getValue<double>("keeplowstate", 0.9998);
        double const keephighstate = arginfo.getValue<double>("keeplowstate", 0.995);
        double const startlowprob = arginfo.getValue<double>("startlowprob", 0.7);

        std::map<int,double> state_erate;
        std::map<int,double> state_erate_sq;
       
        state_erate[state_error_low] = arginfo.getValue<double>("eratelow", 0.15);
        state_erate[state_error_high] = arginfo.getValue<double>("eratehigh",0.25);
        state_erate_sq[state_error_low] = arginfo.getValue<double>("eratelowstddev", 0.03);
        state_erate_sq[state_error_high] = arginfo.getValue<double>("eratehighstddev", 0.04);
    
        libmaus2::fastx::StreamFastAReaderWrapper SFAR(std::cin);
        libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
        std::ostringstream seqstr;
        
        while ( SFAR.getNextPatternUnlocked(pattern) )
            seqstr << pattern.spattern;
        
        std::string const seq = seqstr.str();
        
        std::map< int, std::map<int,double> > state_change_map;
        std::map< int, double > state_start_map;
        
        state_change_map[state_error_low][state_error_low ] = keeplowstate;
        state_change_map[state_error_low][state_error_high] = 1-state_change_map[state_error_low][state_error_low ];
        state_change_map[state_error_high][state_error_high] = keephighstate;
        state_change_map[state_error_high][state_error_low] = 1-state_change_map[state_error_high][state_error_high ];
        state_start_map[state_error_low] = startlowprob;
        state_start_map[state_error_high] = 1 - state_start_map[state_error_low];
        
        uint64_t runid = 0;
        uint64_t readid = 0;

        for ( uint64_t trav = 0; trav < numtraversals; ++trav )
        {
            uint64_t pp = 0;
            typedef std::pair<uint64_t,uint64_t> upair;
            std::vector<upair> poslenvec;
            
            bool const strand = libmaus2::random::Random::rand8() & 1;
            
            while ( pp < seq.size() )
            {
                int64_t len = libmaus2::random::GaussianRandom::random(readlenstddev,readlenavg);
                
                len = std::min(len,static_cast<int64_t>(seq.size()-pp));
                
                poslenvec.push_back(upair(pp,len));
                
                pp += len;
            }
        
            for ( uint64_t z = 0; z < poslenvec.size(); ++z )
            {
                if ( libmaus2::random::UniformUnitRandom::uniformUnitRandom() < droprate )
                    continue;
            
                std::ostringstream baseostr;
                std::ostringstream errostr;
                
                uint64_t const len = poslenvec[z].second;
                uint64_t const pos = strand ? poslenvec[z].first : (seq.size() - poslenvec[z].first - len);

                std::string sub = strand ? seq.substr(pos,len) : libmaus2::fastx::reverseComplementUnmapped(seq.substr(pos,len));
                
                std::vector < state_enum > states;
                state_enum state = (libmaus2::random::UniformUnitRandom::uniformUnitRandom() < state_start_map[state_error_low]) ? state_error_low : state_error_high;
                for ( uint64_t i = 0; i < sub.size(); ++i )
                {
                    states.push_back(state);
                    
                    if ( libmaus2::random::UniformUnitRandom::uniformUnitRandom() < state_change_map[state][state_error_low] )
                        state = state_error_low;
                    else
                        state = state_error_high;
                }
                
                std::map<uint64_t,std::vector<err_type_enum> > errM;
                for ( uint64_t i = 0; i < sub.size(); ++i )
                    errM[i] = std::vector<err_type_enum>(0);
                
                
                uint64_t low = 0;
                while ( low < sub.size() )
                {
                    uint64_t high = low;
                    while ( high < sub.size() && states[high] == states[low] )
                        ++high;

                    double erate = -1;
                    while ( (erate = libmaus2::random::GaussianRandom::random(state_erate_sq[states[low]],state_erate[states[low]])) < 0 )
                    {
                    
                    }
                    
                    uint64_t const numerr = std::floor(erate * (high-low) + 0.5);
                        
                    errostr << "einterval(" << "[" << low << "," << high << ")" << ",erate=" << erate << ",numerr=" << numerr << ')' << ';';
                    
                    uint64_t errplaced = 0;
                    
                    while ( errplaced < numerr )
                    {
                        double const p = libmaus2::random::UniformUnitRandom::uniformUnitRandom();
                        uint64_t const errpos = std::floor(libmaus2::random::UniformUnitRandom::uniformUnitRandom() * (high-low-1) + 0.5) + low;
                        
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
                    
                    low = high;
                }
             
                errostr << "pos=" << pos << ';';
                errostr << "strand=" << strand << ';';
                
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
                            errostr << 'i' << static_cast<char>(insbase);
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
                            errostr << 's' << static_cast<char>(insbase);
                        }
                        else
                        {
                           baseostr.put(sub[pos]); 
                           errostr << 'o';
                        }
                    }
                    else
                    {
                        errostr.put('d');
                    }
                }
                
                std::cout << '>' << 'L' << runid << '/' << (readid++) << '/' << 0 << '_' << baseostr.str().size() << " RQ=0.851 " << errostr.str() << '\n';

                uint64_t b_low = 0;
                std::string const bases = baseostr.str();
                
                while ( b_low < bases.size() )
                {
                    uint64_t const high = std::min(b_low + 80, static_cast<uint64_t>(bases.size()));
                    
                    std::cout.write(bases.c_str() + b_low, high-b_low);
                    std::cout.put('\n');
                
                    b_low = high;
                }
            }
        }
        
    }
    catch(std::exception const & ex)
    {
        std::cerr << ex.what() << std::endl;
        return EXIT_FAILURE;
    }
}
