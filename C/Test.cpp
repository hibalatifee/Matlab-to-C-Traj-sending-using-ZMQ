#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <zmq.h>
#include <iostream>
#include <fstream>

using namespace std;



int main (int argc, char *argv[])
{
		double px=0;
		double py=0;
		double pz=0;
		char message [30001]={};
		void *context1 = zmq_ctx_new();
	

    //  Connect our subscriber socket
	void *subscriber = zmq_socket(context1, ZMQ_SUB);
	
	//sprintf_s(REQUEST,"%0.3f", 1.5);
    //zmq_send (subscriber, REQUEST, 30, 0);
    //zmq_connect(subscriber, "tcp://localhost:5566");

    //  Synchronize with publisher
    //void *sync = zmq_socket (context1, ZMQ_PUSH);
    //zmq_connect(sync, "tcp://localhost:5567");
    //zmq_send (sync, "", 30, 0);
	cout <<"Connecting to the publisher"<< endl;
	zmq_connect(subscriber, "tcp://localhost:5577");
	zmq_setsockopt(subscriber, ZMQ_SUBSCRIBE,"",0);

    //  Get updates, expect random Ctrl-C death
    while (1) {
		
		cout <<"Message:"<< endl;
		zmq_recv (subscriber, message, 3001, 0);
		cout << message<< endl;

		char *token;
		char *rest=message;

	token=strtok_s(rest, " ", &rest);
	px=atof(token);

	token=strtok_s(rest, " ", &rest);
	py=atof(token);
	    
	token=strtok_s(rest, " ", &rest);
	pz=atof(token);

    cout<<"Recieved Position Px="<<px<<"Recieved Position Py="<<py<<"Recieved Position Pz="<<pz<<  endl;
        
    }
		zmq_disconnect(subscriber,"tcp://localhost:5575");
		zmq_close (subscriber);
		//zmq_close (sync);
		zmq_ctx_destroy (context1);
		zmq_ctx_term (context1);
    return 0;
}