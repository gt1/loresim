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
#include <libmaus2/fastx/StreamFastAReader.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
	
		libmaus2::fastx::StreamFastAReaderWrapper in(std::cin);	
		libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
		
		uint64_t const runid = arginfo.getValueUnsignedNumeric<uint64_t>("run",0);
		uint64_t readid = arginfo.getValueUnsignedNumeric<uint64_t>("readidbase",0);
		uint64_t const cols = arginfo.getValueUnsignedNumeric<uint64_t>("cols",80);
		
		while ( in.getNextPatternUnlocked(pattern) )
		{
			std::string const s = pattern.spattern;
			char const * c = s.c_str();

			std::cout << '>' << 'L' << runid << '/' << (readid++) << '/' << 0 << '_' << s.size() << " RQ=0.851\n";
			uint64_t low = 0;
			
			while ( low < s.size() )
			{
				uint64_t const high = std::min(low+cols,static_cast<uint64_t>(s.size()));
				
				std::cout.write(c+low,high-low);
				std::cout.put('\n');
				
				low = high;
			}
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
