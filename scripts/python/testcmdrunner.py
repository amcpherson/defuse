import cmdrun

cmdrun = cmdrun.cmdrun("test","log","sge")
cmdrun.mem = 1
cmdrun.max_parallel = 2
cmdrun.file_timeout = 10

cmdrun.padd(["echo asdf1 | true"], {}, {})
cmdrun.padd(["echo asdf2 | false"], {}, {})
cmdrun.padd(["echo asdf3 | true"], {}, {})
cmdrun.padd(["true"], {}, {})
cmdrun.padd(["touch %(in1)s"], {}, {'in1':"test1"})
cmdrun.padd(["true"], {}, {'in2':"test2"})
#cmdrun.padd("./testcmdrunner2.pl testinner", {}, {})
cmdrun.prun()

