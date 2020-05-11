-module('mergesort').
-import(lists,[map/2,flatten/1,merge/2,split/2]).
-export([main/1,merge_sort/1,wrt/2,parallel_merge_sort/4,rcvarray/1]).

% Note that below implementation is strictly for a fixed number of processes
% In our case 8, hence the parallelisation stops at depth 3

merge_sort(Ar) ->
    N = length(Ar),
    if
        N == 1 ->
            Ar;
        true ->
            {Left_Ar,Right_Ar} = split(N div 2, Ar),
            Sorted_left = merge_sort(Left_Ar),
            Sorted_right = merge_sort(Right_Ar),
            lists:merge(Sorted_left,Sorted_right)
    end.

rcvarray(PID) -> 
    receive
        {PID, Ar} -> Ar
    end.

parallel_merge_sort(Rank,Ar,Depth,Parent_ID) ->
    N = length(Ar),
    % io:format("Rank : ~w\n",[Rank]),
    % io:format("Rank : ~w ~w\n",[Rank,Ar]),
    if
        N < 9 ->
            Parent_ID ! {self(),merge_sort(Ar)};
            % io:format("Rank : ~w work done !!\n",[Rank]);
        true ->
            if
                Depth > 2 ->
                    Parent_ID ! {self(),merge_sort(Ar)};
                    % io:format("Rank : ~w work done !!\n",[Rank]);
                true ->
                    {Left_Ar,Right_Ar} = split(N div 2, Ar),
                    Pid1 = spawn('mergesort', parallel_merge_sort, [Rank*2+1, Left_Ar, Depth+1, self()]),
                    Pid2 = spawn('mergesort', parallel_merge_sort, [Rank*2+2, Right_Ar, Depth+1, self()]),
                    Sorted_left = rcvarray(Pid1),
                    Sorted_right = rcvarray(Pid2),
                    Parent_ID ! {self(), merge(Sorted_left,Sorted_right)}
                    % io:format("Rank : ~w work done !!\n",[Rank])
            end
    end.

wrt(Output,Ar) ->
    N = length(Ar),
    {Rem_Ar, Last} = lists:split(N - 1, Ar),
    {ok,File} = file:open(Output,[write]),
    Text = lists:flatten([io_lib:format("~p ",[Element]) || Element<-Rem_Ar]),
    io:format(File,"~s",[Text]),
    Terminator = io_lib:format("~w~n",Last),
    file:write(File,Terminator),
    file:close(Output).

main(Args) ->

    [Input,Output] = Args,

    {ok,File} = file:read_file(Input),
    Text = string:tokens(erlang:binary_to_list(File), " \n"),
    Ar = lists:map(fun(X) -> list_to_integer(X) end,Text),
    PID = spawn('mergesort',parallel_merge_sort,[0,Ar,1,self()]),
    Sorted_ar = rcvarray(PID),
    wrt(Output,Sorted_ar).
