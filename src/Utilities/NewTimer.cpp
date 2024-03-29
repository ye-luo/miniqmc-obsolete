////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Ken Esler, kpesler@gmail.com,
//    University of Illinois at Urbana-Champaign
// Jeongnim Kim, jeongnim.kim@gmail.com,
//    University of Illinois at Urbana-Champaign
// Jeremy McMinnis, jmcminis@gmail.com,
//    University of Illinois at Urbana-Champaign
// Mark Dewing, mdewing@anl.gov,
//    Argonne National Laboratory
//
// File created by:
// Ken Esler, kpesler@gmail.com,
//    University of Illinois at Urbana-Champaign
////////////////////////////////////////////////////////////////////////////////

/** @file NewTimer.cpp
 * @brief Implements TimerManager
 */
#include "Configuration.h"
#include "Utilities/NewTimer.h"
#include <map>
#include <limits>
#include <cstdio>

namespace qmcplusplus
{
TimerManagerClass TimerManager;

bool timer_max_level_exceeded = false;

void TimerManagerClass::addTimer(NewTimer *t)
{
#pragma omp critical
  {
    if (t->get_name().find(TIMER_STACK_SEPARATOR) != std::string::npos)
    {
      app_log() << "Warning: Timer name (" << t->get_name()
                << ") should not contain the character "
                << TIMER_STACK_SEPARATOR << std::endl;
    }

    if (timer_name_to_id.find(t->get_name()) == timer_name_to_id.end())
    {
      t->set_id(max_timer_id);
      timer_id_name[t->get_id()]      = t->get_name();
      timer_name_to_id[t->get_name()] = t->get_id();
      if (max_timer_id >= std::numeric_limits<timer_id_t>::max())
      {
        max_timers_exceeded = true;
        app_log() << "Number of timers exceeds limit ("
                  << static_cast<int>(std::numeric_limits<timer_id_t>::max())
                  << ").   Adjust timer_id_t in NewTimer.h and recompile."
                  << std::endl;
      }
      else
      {
        max_timer_id++;
      }
    }
    else
    {
      t->set_id(timer_name_to_id[t->get_name()]);
    }
    t->set_manager(this);
    t->set_active_by_timer_threshold(timer_threshold);
    TimerList.push_back(t);
  }
}

NewTimer *TimerManagerClass::createTimer(const std::string &myname,
                                         timer_levels mytimer)
{
  NewTimer *t = new NewTimer(myname, mytimer);
  addTimer(t);
  return t;
}

void TimerManagerClass::reset()
{
  for (int i = 0; i < TimerList.size(); i++)
    TimerList[i]->reset();
}

void TimerManagerClass::set_timer_threshold(const timer_levels threshold)
{
  timer_threshold = threshold;
  for (int i = 0; i < TimerList.size(); i++)
  {
    TimerList[i]->set_active_by_timer_threshold(timer_threshold);
  }
}

void TimerManagerClass::collate_flat_profile(FlatProfileData &p)
{
  for (int i = 0; i < TimerList.size(); ++i)
  {
    NewTimer &timer = *TimerList[i];
    nameList_t::iterator it(p.nameList.find(timer.get_name()));
    if (it == p.nameList.end())
    {
      int ind                      = p.nameList.size();
      p.nameList[timer.get_name()] = ind;
      p.timeList.push_back(timer.get_total());
      p.callList.push_back(timer.get_num_calls());
    }
    else
    {
      int ind = (*it).second;
      p.timeList[ind] += timer.get_total();
      p.callList[ind] += timer.get_num_calls();
    }
  }
}

struct ProfileData
{
  double time;
  double calls;

  ProfileData &operator+=(const ProfileData &pd)
  {
    time += pd.time;
    calls += pd.calls;
    return *this;
  }
};

int get_level(const std::string &stack_name)
{
  int level = 0;
  for (int i = 0; i < stack_name.length(); i++)
  {
    if (stack_name[i] == TIMER_STACK_SEPARATOR)
    {
      level++;
    }
  }
  return level;
}

std::string get_leaf_name(const std::string &stack_name)
{
  int pos = stack_name.find_last_of(TIMER_STACK_SEPARATOR);
  if (pos == std::string::npos)
  {
    return stack_name;
  }

  return stack_name.substr(pos + 1, stack_name.length() - pos);
}

void TimerManagerClass::get_stack_name_from_id(const StackKey &key,
                                               std::string &stack_name)
{
  for (int i = 0; i < StackKey::max_level; i++)
  {
    std::string &timer_name = timer_id_name[key.get_id(i)];
    if (key.get_id(i) == 0) break;
    if (i > 0)
    {
      stack_name += TIMER_STACK_SEPARATOR;
    }
    stack_name += timer_name;
  }
}

void TimerManagerClass::collate_stack_profile(StackProfileData &p)
{
#ifdef USE_STACK_TIMERS
  // Put stacks from all timers into one data structure
  // By naming the timer stacks as 'timer1/timer2', etc, the ordering done by
  // the
  // map's keys will also place the stacks in depth-first order.
  // The order in which sibling timers are encountered in the code is not
  // preserved. They will be ordered alphabetically instead.
  std::map<std::string, ProfileData> all_stacks;
  for (int i = 0; i < TimerList.size(); ++i)
  {
    NewTimer &timer = *TimerList[i];
    std::map<StackKey, double>::iterator stack_id_it =
        timer.get_per_stack_total_time().begin();
    for (; stack_id_it != timer.get_per_stack_total_time().end(); stack_id_it++)
    {
      ProfileData pd;
      const StackKey &key = stack_id_it->first;
      std::string stack_name;
      get_stack_name_from_id(key, stack_name);
      pd.time  = timer.get_total(key);
      pd.calls = timer.get_num_calls(key);

      all_stacks[stack_name] += pd;
    }
  }

  // Fill in the output data structure (but don't compute exclusive time yet)
  std::map<std::string, ProfileData>::iterator si = all_stacks.begin();
  int idx = 0;
  for (; si != all_stacks.end(); ++si)
  {
    std::string stack_name = si->first;
    p.nameList[stack_name] = idx;
    p.names.push_back(stack_name);
    p.timeList.push_back(si->second.time);
    p.timeExclList.push_back(si->second.time);
    p.callList.push_back(si->second.calls);
    idx++;
  }

  // Subtract times of immediate children to get exclusive time
  for (idx = 0; idx < p.timeList.size(); idx++)
  {
    int start_level = get_level(p.names[idx]);
    for (int i = idx + 1; i < p.timeList.size(); i++)
    {
      int level = get_level(p.names[i]);
      if (level == start_level + 1)
      {
        p.timeExclList[idx] -= p.timeExclList[i];
      }
      if (level == start_level)
      {
        break;
      }
    }
  }
#endif
}

void TimerManagerClass::print()
{
#if ENABLE_TIMERS
#ifdef USE_STACK_TIMERS
  printf("Stack timer profile\n");
  print_stack();
#else
  printf("\nFlat profile\n");
  print_flat();
#endif
#endif
}

void TimerManagerClass::print_flat()
{
#if ENABLE_TIMERS
  FlatProfileData p;

  collate_flat_profile(p);

  {
#pragma omp master
    {
      std::map<std::string, int>::iterator it(p.nameList.begin()),
          it_end(p.nameList.end());
      while (it != it_end)
      {
        int i = (*it).second;
        // if(callList[i]) //skip zeros
        printf("%-40s  %9.4f  %13ld  %16.9f  %12.6f TIMER\n",
               (*it).first.c_str(), p.timeList[i], p.callList[i],
               p.timeList[i] / (static_cast<double>(p.callList[i]) +
                                std::numeric_limits<double>::epsilon()),
               p.timeList[i] / static_cast<double>(omp_get_max_threads()));
        ++it;
      }
    }
  }
#endif
}

void pad_string(const std::string &in, std::string &out, int field_len)
{
  int len     = in.size();
  int pad_len = std::max(field_len - len, 0);
  std::string pad_str(pad_len, ' ');
  out = in + pad_str;
}

void TimerManagerClass::print_stack()
{
#if ENABLE_TIMERS
  StackProfileData p;

  collate_stack_profile(p);

  {
    if (timer_max_level_exceeded)
    {
      printf("Warning: Maximum stack level (%d) exceeded.  Results may be "
             "incorrect.\n",
             StackKey::max_level);
      printf("Adjust StackKey in NewTimer.h and recompile.\n");
    }

    int indent_len   = 2;
    int max_name_len = 0;
    for (int i = 0; i < p.names.size(); i++)
    {
      std::string stack_name = p.names[i];
      int level              = get_level(stack_name);
      std::string name       = get_leaf_name(stack_name);
      int name_len           = name.size() + indent_len * level;
      max_name_len           = std::max(name_len, max_name_len);
    }

    std::string timer_name;
    pad_string("Timer", timer_name, max_name_len);
    printf("%s  %-9s  %-9s  %-10s  %-13s\n", timer_name.c_str(),
           "Inclusive_time", "Exclusive_time", "Calls", "Time_per_call");
    for (int i = 0; i < p.names.size(); i++)
    {
      std::string stack_name = p.names[i];
      int level              = get_level(stack_name);
      std::string name       = get_leaf_name(stack_name);
      std::string indent_str(indent_len * level, ' ');
      std::string indented_str = indent_str + name;
      std::string padded_name_str;
      pad_string(indented_str, padded_name_str, max_name_len);
      printf("%s  %9.4f  %9.4f  %13ld  %16.9f\n", padded_name_str.c_str(),
             p.timeList[i], p.timeExclList[i], p.callList[i],
             p.timeList[i] / (static_cast<double>(p.callList[i]) +
                              std::numeric_limits<double>::epsilon()));
    }
  }
#endif
}

// Might want some sort of structured output for timing data - either xml or
// yaml
#if 0
void
TimerManagerClass::output_timing(Communicate *comm, Libxml2Document &doc, xmlNodePtr root)
{
#if ENABLE_TIMERS
#ifdef USE_STACK_TIMERS
  StackProfileData p;

  collate_stack_profile(comm, p);

  if(comm == NULL || comm->rank() == 0)
  {
    xmlNodePtr timing_root = doc.addChild(root, "timing");
    doc.addChild(timing_root, "max_stack_level_exceeded", timer_max_level_exceeded?"yes":"no");
    doc.addChild(timing_root, "max_timers_exceeded", max_timers_exceeded?"yes":"no");
    std::vector<xmlNodePtr> node_stack;
    node_stack.push_back(timing_root);
    xmlNodePtr current_root = timing_root;

    for (int i = 0; i < p.names.size(); i++)
    {
      std::string stack_name = p.names[i];
      int level = get_level(stack_name);
      std::string name = get_leaf_name(stack_name);

      std::string indent_str(2*level, ' ');

      xmlNodePtr timer = doc.addChild(current_root, "timer");
      doc.addChild(timer, "name", name);
      doc.addChild(timer, "time_incl", p.timeList[i]);
      doc.addChild(timer, "time_excl", p.timeExclList[i]);
      doc.addChild(timer, "calls", p.callList[i]);

      int next_level = level;
      if (i+1 < p.names.size())
      {
        next_level = get_level(p.names[i+1]);
      }

      if (next_level > level)
      {
        xmlNodePtr next_node = doc.addChild(timer, "includes");
        node_stack.push_back(next_node);
        current_root = next_node;
      }
      if (next_level < level)
      {
        for (int j = 0; j < level-next_level; j++)
        {
          node_stack.pop_back();
          current_root = node_stack.back();
        }
      }
    }
  }

#endif
#endif
}
#endif

void NewTimer::set_active_by_timer_threshold(const timer_levels threshold)
{
  if (timer_level <= threshold)
  {
    active = true;
  }
  else
  {
    active = false;
  }
}
}
