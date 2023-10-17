library IEEE;
use IEEE.STD_LOGIC_1164.ALL;

entity snake_tb is
--  Port ( );
end snake_tb;

architecture Behavioral of snake_tb is

    component snake is
        Port ( clk : in STD_LOGIC;
               rst_n : in STD_LOGIC;
               dir : in STD_LOGIC;
               en : in STD_LOGIC;
               led : out STD_LOGIC_VECTOR (15 downto 0));
    end component;
    
    signal clk: std_logic := '0';
    signal rst_n: std_logic := '1';
    signal dir: std_logic := '0';
    signal en: std_logic := '0';
    signal led: std_logic_vector(15 downto 0);

    begin
        -- Instantiate the component
        snake_i: snake port map(clk => clk, rst_n => rst_n, led => led, dir => dir, en => en);
        -- Define clock behvaiour
        clk <= not clk after 5ns;
        -- Test process
        process begin
            -- Wait for clock for synchronization
            wait for 5ns;
            
            rst_n <= '0'; en <= '0'; dir <= '0';
            wait for 10ns;
            
            rst_n <= '1'; en <= '0'; dir <= '0';
            wait for 50ns;
            
            rst_n <= '1'; en <= '1'; dir <= '0';
            wait for 150ns;
            
            rst_n <= '1'; en <= '1'; dir <= '1';
            wait for 150ns;
            
            rst_n <= '1'; en <= '0'; dir <= '1';
            wait;
        end process;


end Behavioral;
