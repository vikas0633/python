#12t.py - script for making mySQl-add column batch file
import sys

table_name=sys.argv[1]
data=sys.argv[2]
annotation=((sys.argv[3]).split('/'))[-1]


table_name=table_name.replace('normalized_score_cutoff','score') 
print str("DROP TABLE IF EXISTS temp;")
print str("DROP TABLE IF EXISTS test;")
print str("CREATE TABLE temp ("+"`Sequence_"+annotation+"` varchar (200) ,"+" `"+annotation+"` VARCHAR(1000) );")
print str("LOAD DATA LOCAL INFILE '"+ str(data)+ "' INTO TABLE `"+"temp"+"` (`Sequence_"+annotation+"`"+", `"+annotation+"`);")
print str("CREATE TABLE test SELECT * FROM `"+((table_name).split('/'))[-1]+"`, temp Where `"+((table_name).split('/'))[-1]+"`.`#Sequence`=temp.`Sequence_"+annotation+"`;")
print str("drop table `"+((table_name).split('/'))[-1]+"`;")
print str("RENAME TABLE test TO `"+((table_name).split('/'))[-1]+"`;")
print str("ALTER TABLE `"+ ((table_name).split('/'))[-1] +"` DROP " +"`Sequence_"+annotation+"`;" )
print str("DROP TABLE temp;")