library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity task1_tb is
    --  Port ( );
end task1_tb;

architecture Behavioral of task1_tb is

    component task1 is
        Port ( x : in STD_LOGIC;
               y : in STD_LOGIC;
               z : in STD_LOGIC;
               f1 : out STD_LOGIC;
               f2 : out STD_LOGIC);
    end component;
    
    signal x: std_logic := '0';
    signal y: std_logic := '0';
    signal z: std_logic := '0';
    signal f1: std_logic;
    signal f2: std_logic;
    
    begin
        task1_i: task1
        port map(x => x, y => y, z => z, f1 => f1, f2 => f2);
        
        -- Go through the Karnaugh map value by value
        process begin
            -- 0
--          x <= '0';
--          y <= '0';
--          z <= '0';
            wait for 50ns;
            -- 1
            y <= '1';
            wait for 50ns;
            -- 3
            z <= '1';
            wait for 50ns;
            -- 2
            y <= '0';
            wait for 50ns;
            -- 4
            x <= '1';
            z <= '0';
            wait for 50ns;
            -- 5
            y <= '1';
            wait for 50ns;
            -- 6
            z <= '1';
            wait for 50ns;
            -- 7
            y <= '0';
            wait;
        end process;

end Behavioral;
