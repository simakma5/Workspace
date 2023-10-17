library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity task2_tb is
--  Port ( );
end task2_tb;

architecture Behavioral of task2_tb is

    component task2 is
        Port ( D : in STD_LOGIC;
               C : in STD_LOGIC;
               B : in STD_LOGIC;
               A : in STD_LOGIC;
               f : out STD_LOGIC);
    end component;
    
    signal D: std_logic := '0';
    signal C: std_logic := '0';
    signal B: std_logic := '0';
    signal A: std_logic := '0';
    signal f: std_logic;

    begin
        task2_i: task2 port map(D => D, C => C, B => B, A => A, f => f);
        
         -- Go through the Karnaugh map value by value
        process begin
            -- 0
--          D <= '0';
--          C <= '0';
--          B <= '0';
--          A <= '0';
            wait for 50ns;
            -- 1
            A <= '1';
            wait for 50ns;
            -- 3
            B <= '1';
            wait for 50ns;
            -- 2
            A <= '0';
            wait for 50ns;
            -- 4
            C <= '1';
            B <= '0';
            wait for 50ns;
            -- 5
            A <= '1';
            wait for 50ns;
            -- 7
            B <= '1';
            wait for 50ns;
            -- 6
            A <= '0';
            wait for 50ns;
            -- 12
            D <= '1';
            B <= '0';
            wait for 50ns;
            -- 13
            A <= '1';
            wait for 50ns;
            -- 15
            B <= '1';
            wait for 50ns;
            -- 14
            A <= '0';
            wait for 50ns;
            -- 8
            C <= '0';
            B <= '0';
            wait for 50ns;
            -- 9
            A <= '1';
            wait for 50ns;
            -- 11
            B <= '1';
            wait for 50ns;
            -- 12
            A <= '0';
            wait;
        end process;


end Behavioral;
