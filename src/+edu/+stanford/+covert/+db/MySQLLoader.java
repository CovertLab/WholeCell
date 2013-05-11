//package edu.stanford.covert.db;

import com.mysql.jdbc.Driver;
import java.sql.Connection;
import java.sql.DriverManager;

public class MySQLLoader
{
    public static Connection makeConnection(String hostName, String schema, String userName, String password)
        throws Exception
    {
        DriverManager.registerDriver(new com.mysql.jdbc.Driver());
        return DriverManager.getConnection("jdbc:mysql://"+hostName+"/"+schema, userName, password);
    }
}