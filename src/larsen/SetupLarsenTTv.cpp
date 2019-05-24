#include "SetupLarsenTTv.h"
#include "DataLarsenTTv.h"
#include "LarsenTTv.h"

// Old LARSEN logo
//
// SetupLarsenTTv::SetupLarsenTTv(){
//       std::cout << "                                                                    " << std::endl;
//       std::cout << "              _ _ __ __ ________   __         __ \\  |  /           " << std::endl;
//       std::cout << "     _ _ __ __  _______________   /  \\_/\\ /\\_/  \\       -       " << std::endl;
//       std::cout << "              _ _ __ __ ________  \\__/   V   \\__/     __          " << std::endl;
//       std::cout << "                                                 /   /  \\          " << std::endl;
//       std::cout << "   _      _   ___  ___ ___ _  _                      \\__/          " << std::endl;
//       std::cout << "  | |    /_\\ | _ \\/ __| __| \\| |  LAgrangian               \\    " << std::endl;
//       std::cout << "  | |__ / _ \\|   /\\__ \\ _|| .` |  Reactor for StrEams  \\  \\ \\ " << std::endl;
//       std::cout << "  |____/_/ \\_\\_|_\\|___/___|_|\\_|  in Nonequilibrium     \\  \\  " << std::endl;
//       std::cout << "                                                            \\      " << std::endl;
//       std::cout << "                                                                    " << std::endl;
// }


SetupLarsenTTv::SetupLarsenTTv(){
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

Data* SetupLarsenTTv::getData(Mutation::Mixture& l_mix,
                             const std::vector< std::string > l_input_file)
{ 
  return new DataLarsenTTv(l_mix, l_input_file); 
}

Problem* SetupLarsenTTv::getProblem(Mutation::Mixture& l_mix, Data& l_data)
{ 
  return new LarsenTTv(l_mix, l_data); 
}

