/*
 * All of the documentation and software included in the
 * Alchemy Software is copyrighted by Stanley Kok, Parag
 * Singla, Matthew Richardson, Pedro Domingos, Marc
 * Sumner, Hoifung Poon, Daniel Lowd, and Jue Wang.
 * 
 * Copyright [2004-09] Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, Daniel Lowd, and Jue Wang. All rights reserved.
 * 
 * Contact: Pedro Domingos, University of Washington
 * (pedrod@cs.washington.edu).
 * 
 * Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that
 * the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the
 * following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the
 * above copyright notice, this list of conditions and the
 * following disclaimer in the documentation and/or other
 * materials provided with the distribution.
 * 
 * 3. All advertising materials mentioning features or use
 * of this software must display the following
 * acknowledgment: "This product includes software
 * developed by Stanley Kok, Parag Singla, Matthew
 * Richardson, Pedro Domingos, Marc Sumner, Hoifung
 * Poon, Daniel Lowd, and Jue Wang in the Department of
 * Computer Science and Engineering at the University of
 * Washington".
 * 
 * 4. Your publications acknowledge the use or
 * contribution made by the Software to your research
 * using the following citation(s): 
 * Stanley Kok, Parag Singla, Matthew Richardson and
 * Pedro Domingos (2005). "The Alchemy System for
 * Statistical Relational AI", Technical Report,
 * Department of Computer Science and Engineering,
 * University of Washington, Seattle, WA.
 * http://alchemy.cs.washington.edu.
 * 
 * 5. Neither the name of the University of Washington nor
 * the names of its contributors may be used to endorse or
 * promote products derived from this software without
 * specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE UNIVERSITY OF WASHINGTON
 * AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE UNIVERSITY
 * OF WASHINGTON OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
#ifndef STRINT_H_JUN_27_2005
#define STRINT_H_JUN_21_2005

#include <cstdlib>
#include <ext/hash_map>
//using namespace __gnu_cxx;


struct StrInt
{
 public:
  StrInt() : str_(NULL), ints_(NULL) {}

  StrInt(const char* const & s, const int& i) 
  {
    str_ = new char[strlen(s)+1];
    strcpy(str_,s);
    ints_ = new Array<int>;
    ints_->append(i);
  }

  StrInt(const char* const & s) 
  {
    str_ = new char[strlen(s)+1];
    strcpy(str_,s);
    ints_ = NULL;
  }

  ~StrInt() { delete [] str_; delete ints_; }

  void addInt(const int& i)
  {
    ints_->append(i);
  }

  char* str_;
  Array<int>* ints_;
};


class HashStrInt
{
 public:
  size_t operator()(const StrInt* const & s) const
  {
    return __gnu_cxx::hash<char const *>()(s->str_);
  }
};


class EqualStrInt
{
 public:
  bool operator()(const StrInt* const & s1, const StrInt* const & s2) const
  {
    return strcmp(s1->str_, s2->str_) == 0;
  }
};



#endif
