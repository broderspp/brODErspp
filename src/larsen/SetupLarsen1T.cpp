#include "SetupLarsen1T.h"
#include "DataLarsen1T.h"
#include "Larsen1T.h"

// Old logo
//SetupLarsen1T::SetupLarsen1T(){
//      std::cout << "                                                                    " << std::endl;
//      std::cout << "              _ _ __ __ ________   __         __ \\  |  /           " << std::endl;
//      std::cout << "     _ _ __ __  _______________   /  \\_/\\ /\\_/  \\       -       " << std::endl;
//      std::cout << "              _ _ __ __ ________  \\__/   V   \\__/     __          " << std::endl;
//      std::cout << "                                                 /   /  \\          " << std::endl;
//      std::cout << "   _      _   ___  ___ ___ _  _                      \\__/          " << std::endl;
//      std::cout << "  | |    /_\\ | _ \\/ __| __| \\| |  LAgrangian               \\    " << std::endl;
//      std::cout << "  | |__ / _ \\|   /\\__ \\ _|| .` |  Reactor for StrEams  \\  \\ \\ " << std::endl;
//      std::cout << "  |____/_/ \\_\\_|_\\|___/___|_|\\_|  in Nonequilibrium     \\  \\  " << std::endl;
//      std::cout << "                                                            \\      " << std::endl;
//      std::cout << "                                                                    " << std::endl;
//}

SetupLarsen1T::SetupLarsen1T(){
  std::cout << "                                                                                                     " << std::endl;
  std::cout << "      /\\                                _______ ___                                                 " << std::endl; 
  std::cout << "     /  \\         ____________          \\      \\  /                                               " << std::endl; 
  std::cout << "    /    \\       /            \\____      \\      \\/0_______        __   ___   ___  _________  __  " << std::endl;      
  std::cout << "    V\\  /V      /   .........      \\     (      0,`/ / /  \\      / /  / _ | / _ \\/ __/ __/ |/ /  " << std::endl;     
  std::cout << "      \\ \\     __\\    .       ....   \\    (                o)    / /__/ __ |/ , _/\\ \\/ _//    / " << std::endl;   
  std::cout << "       ) \\   /  \\)   .  /\\___   ..   \\__/ .                |   /____/_/ |_/_/|_/___/___/_/|_/    " << std::endl;    
  std::cout << "       (  \\ / ....  .  (     )  ..........       /wwwwwwwww'                                        " << std::endl;     
  std::cout << "        \\  v .  _ ..    )    (   ..........     (  ---------------------------------------------->  " << std::endl;    
  std::cout << "         ). .  / \\_____/      \\   ..  _____      \\MMMMMMMM\\   LAgrangian Reactor for StrEams     " << std::endl;  
  std::cout << "         \\ .  /                \\_____/     \\______________/                     in Nonequilibrium " << std::endl; 
  std::cout << "          \\__/                                                                                      " << std::endl; 
  std::cout << "                                                                                                     " << std::endl;
}

Data* SetupLarsen1T::getData(Mutation::Mixture& l_mix, 
                             const std::vector< std::string > l_input_file)
{ 
  return new DataLarsen1T(l_mix, l_input_file); 
}

Problem* SetupLarsen1T::getProblem(Mutation::Mixture& l_mix, Data& l_data)
{ 
  return new Larsen1T(l_mix, l_data); 
}

