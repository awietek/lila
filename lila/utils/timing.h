// Copyright 2018 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LILA_TIMING_H_
#define LILA_TIMING_H_

#include <chrono>

#define LILA_CLK(t) t=std::chrono::high_resolution_clock::now()
#define LILA_TIME(label, t1, t2) std::cout<<label<<" "<<std::chrono::duration_cast<std::chrono::seconds>(t2-t1).count()<<"\n"
#define LILA_TIME_MILLI(label, t1, t2) std::cout<<label<<" "<<std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count()<<"\n"
#define LILA_TIME_MICRO(label, t1, t2) std::cout<<label<<" "<<std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()<<"\n"
#define LILA_TIME_NANO(label, t1, t2) std::cout<<label<<" "<<std::chrono::duration_cast<std::chrono::nanoseconds>(t2-t1).count()<<"\n"

#endif
