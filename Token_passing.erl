-module('Token_passing').
-export([main/1, initiator/3, simulate/5, rcv_wrt/1]).

rcv_wrt(Output) ->
	receive
		{Received_Tuple} ->
			[Recvby,Token,Sentby] = Received_Tuple,
			if
				Sentby == 0 ->
					{ok, File} = file:open(Output,[write]);
				true ->
					{ok, File} = file:open(Output,[append])
			end,
			Text = io_lib:format("Process ~w received token ~w from process ~w.~n", [Recvby,Token,Sentby]),
			file:write(File,Text)
	end.

simulate(Rank, Numprocs, Token, Output, Root_ID) ->
	% SID = self(),
	if
		Rank == 0 ->
			ok;
		Rank > 0 ->
			rcv_wrt(Output)
	end,

	if
		Rank < Numprocs-1 ->
			PID = spawn('Token_passing', simulate, [Rank + 1, Numprocs, Token, Output, Root_ID]),
			% if
			% 	(SID == PID) ->
			% 		PID ! {[Rank+1, Token+1, Rank]};
			% 	true ->
			PID ! {[Rank+1, Token, Rank]};
			% end;
		Rank == Numprocs-1 ->
			Root_ID ! {[0, Token, Rank]}
	end.

initiator(Numprocs,Token,Output) ->
	spawn('Token_passing', simulate, [0, Numprocs, Token, Output, self()]),
	rcv_wrt(Output).

main(Args) -> 
	
	[Input,Output] = Args,
	{ok, File} = file:read_file(Input),
	Text = string:tokens(erlang:binary_to_list(File), " \n"),
	[N,T] = Text,
	Numprocs = list_to_integer(N),
	Token = list_to_integer(T),

	initiator(Numprocs,Token,Output).

	
