//
//  ExecutionService.cpp
//  ParallelBalls
//
//  Created by msmith on 9/9/14.
//  Copyright (c) 2014 paluchlab. All rights reserved.
//

#include "ExecutionService.h"

void ExecutionThread::execution(){
        
    while(running || tasks.size()>0){
        std::unique_lock<std::mutex> lk(m);
        if(tasks.size()>0){
            std::function<void()> cmd = tasks.front();
            
            tasks.pop();
            lk.unlock();
            cmd();
        } else{
            cv.wait(lk);
        }
    }
    std::cout<<"thread finished!\n";
}

void ExecutionThread::start(){
        running=true;
        worker = new std::thread([&](){execution();});
    }
    
void ExecutionThread::submit(std::function<void ()> fun){
        std::lock_guard<std::mutex> lk(m);
        tasks.push(fun);
        cv.notify_one();
    }

ExecutionService::ExecutionService(int threads){
    for(int i = 0; i<threads; i++){
        ExecutionThread* t = new ExecutionThread();
        t->start();
        executors.push_back(t);
    }
    next = 0;
    MAX = threads;
}

void ExecutionService::submit(const std::function<void()> &callable){
    std::promise<int>* promise = new std::promise<int>();
    executors[next]->submit( [callable, promise]{callable(); promise->set_value(0);} );
    
    next = (next + 1)%MAX;
    
    promises.push(promise);
    
}

void ExecutionService::waitForExecution(){
    while(!promises.empty()){
        auto f = promises.front();
        promises.pop();
        f->get_future().get();
        delete f;
    }
}

