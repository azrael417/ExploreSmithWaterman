//include kokkos stuff
#ifndef _TYPES
#define _TYPES

#include <string>
#include <kokkos/core/src/Kokkos_Core.hpp>

template<typename T, typename... Args>
using View1D = Kokkos::View<T*, Args...>;

template<typename T, typename... Args>
using View2D = Kokkos::View<T**, Args...>;


template <typename V>
void ViewToString(std::string& result, V view){
  static_assert( V::rank == 1 , "Error: the view has to be one-dimensional.");
  static_assert( std::is_same<char,typename V::value_type>::value, "Error: the view-type has to be char.");
  result.clear();
  for(uint64_t l=0; l<view.extent(0); l++){
    if(view(l)=='\0') break;
    result+=view(l);
  }
}

template <typename V>
void StringToView(V view, const std::string input){
  static_assert( V::rank == 1 , "Error: the view has to be one-dimensional.");
  static_assert( std::is_same<char,typename V::value_type>::value, "Error: the view-type has to be char.");
  for(uint64_t l=0; l<input.size(); l++){
    view(l) = input[l];
  }
  if(input.size() < view.extent(0)) view(input.size())='\0';
}

// range interation policy
typedef Kokkos::MDRangePolicy<Kokkos::Rank<2,Kokkos::Iterate::Right,Kokkos::Iterate::Right>> t_policy;

#endif
